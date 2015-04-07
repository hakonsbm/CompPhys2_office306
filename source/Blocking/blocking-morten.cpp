#include <mpi.h>
#include <armadillo>
#include <math.h>
#include <iostream>
#include <sys/time.h>
#include <fstream>

class Blocking {
public:
    Blocking(std::string path,
            std::string filename,
            int maxb,
            int minb,
            int n_b,
            int node,
            int n_nodes);

    double estimate_error();

private:

    arma::vec local_block;

    int min_block_size;
    int max_block_size;
    int n_block_samples;

    int n_c;

    bool is_master;
    int node;
    int n_nodes;

    std::string filename;
    std::string path;
    std::ofstream file;

    arma::vec data;

    void block_data(int block_size, double &var, double &mean);
    void get_unique_blocks(arma::Row<int> & block_sizes, int & n);
    void get_initial_error();

    double combine_mean(double mean, int n, int n_tot);
    double combine_variance(double var, double mean = 0, int n = 0);

};

Blocking::Blocking(std::string path, std::string filename,
        int maxb,
        int minb,
        int n_b,
        int node,
        int n_nodes) {

    this->filename = filename;
    this->path = path;

    this->node = node;
    this->n_nodes = n_nodes;

    is_master = (node == 0);

    n_block_samples = n_b;
    max_block_size = maxb / n_nodes;
    min_block_size = minb / n_nodes;
    
    //Load the previously stored data (unique data pr. MPI process)
    data.load(path + filename);
    
    n_c = data.n_elem;

    if (is_master) {
        if (max_block_size > this->n_c / 2) {
            std::cout << "invalid local max block size " << max_block_size << std::endl;
            std::cout << "max block size must not be greater than " << n_nodes * this->n_c / 2 << std::endl;
            exit(1);
        }

        if (min_block_size < 1) {
            std::cout << "invalid local min block size " << min_block_size << std::endl;
            std::cout << "min block size must not be lower than n_nodes=" << n_nodes << std::endl;
            exit(1);
        }

        if (n_block_samples > max_block_size - min_block_size) {
            std::cout << "invalid amount of block samples " << n_block_samples << std::endl;
            std::cout << "block samples must be lower or equal " << max_block_size - min_block_size << std::endl;
            exit(1);
        }
    }

    //The outputfile is opened by the master node.
    if (is_master) this->file.open((path + "blocking_test_out.dat").c_str());

}

double Blocking::estimate_error() {

    int block_size, n;
    double error, var, mean;

    //Should match the one from your original calculation (no blocking)
    //Note: ONLY if you are using 'sampling variance'. If not, it will be
    //somewhat different.
    get_initial_error();

    //get_unique_blocks sheds out any duplicate size due to integer division.
    arma::Row<int> block_sizes = arma::zeros<arma::Row<int> >(n_block_samples);
    get_unique_blocks(block_sizes, n);

    //Loop over block sizes
    for (int j = 0; j < n; j++) {
        block_size = block_sizes(j);

        //local result from each block.
        block_data(block_size, var, mean);
        
        //Fancy mathematics to combine correlated variances.
        var = combine_variance(var, mean);

        //Master calculates roots and dumps stuff to file.
        if (is_master) {
            error = sqrt(var / ((n_nodes * n_c) / block_size - 1.0));

            file << block_size * n_nodes << "\t" << error << std::endl;
            if (j % 9 == 0) {
                std::cout << "\rBlocking progress: " << (double) (j + 1) / n * 100 << "%";
                std::cout.flush();
            }
        }
    }

    if (is_master) std::cout << " Done." << std::endl;
    return error;
}

void Blocking::block_data(int block_size, double &var, double &mean) {
    using namespace arma;

    double block_mean;


    mean = 0;
    double mean2 = 0;

    //n_b = number of blocks
    int n_b = n_c / block_size;
    for (int j = 0; j < n_b; j++) {
        
        //averaging the intervals which split the full interval into n_b blocks.
        block_mean = arma::mean(data(span(j*block_size, (j + 1) * block_size - 1)));

        mean += block_mean;
        mean2 += block_mean*block_mean;

    }

    mean /= n_b;

    //Fancy 'sampling variance' since local block sizes can get small,
    //and hence this is a better estimate.
    var = mean2 / (n_b - 1) - n_b * mean * mean / (n_b - 1);

}

void Blocking::get_initial_error() {

    double var = arma::var(data);
    double mean = arma::mean(data);

    var = combine_variance(var, mean);

    if (is_master) std::cout << "Initial stddev: " << sqrt(var / (n_nodes * n_c - 1)) << std::endl;
}

void Blocking::get_unique_blocks(arma::Row<int>& block_sizes, int& n) {

    int block_step_length = (max_block_size - min_block_size) / (n_block_samples - 1);

    for (int j = 0; j < n_block_samples; j++) {
        block_sizes(j) = min_block_size + j * block_step_length;
    }

    block_sizes = arma::unique(block_sizes);

    n = block_sizes.n_elem;

}

double Blocking::combine_mean(double mean, int n, int n_tot) {

    mean *= n;

    MPI_Allreduce(MPI_IN_PLACE, &mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    mean /= n_tot;

    return mean;

}

double Blocking::combine_variance(double var, double mean, int n) {

    if (n == 0) {
        n = data.n_elem;
    }


    int n_tot;
    MPI_Allreduce(&n, &n_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


    double combined_mean = combine_mean(mean, n, n_tot);
    double combined_var = n * (mean - combined_mean)*(mean - combined_mean) + (n - 1) * var;

    if (is_master) {
        MPI_Reduce(MPI_IN_PLACE, &combined_var, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        var = combined_var / (n_tot - 1);
    } else {
        MPI_Reduce(&combined_var, new double(), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }


    return var;
}

double f(double x){
    return 2*x;
}

void integrate_me(std::string path, std::string name, int node, int n_nodes){
    using namespace arma;
    
    double x, local_I;
    
    //unique seed
    std::srand(-time(NULL)+node);
    
    int n_c = 1000000;
    vec blocking_data = zeros<vec>(n_c);

    double I = 0;
    
    int a = 1;
    int b = 2;
    
    for (int i = 0; i < n_c; i++){
        
        x = a + (b-a)*as_scalar(randu<vec>(1));
 
        local_I = (b-a)*f(x);
  
        //store the "local results"
        blocking_data(i) = local_I;
    
        I += local_I;
        
    }
    
    sleep(node);
    cout << "node " << node << " calculated std: " << sqrt(var(blocking_data)/(n_c-1)) << endl;
    
   
    //Slam the local results to file.
    blocking_data.save(path + name);
    blocking_data.reset();
    
    //Reduce and present
    double I_tot = 0;
    MPI_Reduce(&I, &I_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    I_tot /= (n_nodes*n_c);
    
    if (node==0) std::cout << "result: " << I_tot << std::endl;
    
}

int main() {

    using namespace std;

    int node, n_nodes;

    //Initializing MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &node);

    //Unique name for each MPI process
    stringstream name;
    name << "blocking" << node << ".arma";
    string path = "/home/jorgmeister/scratch/testBlocking/";

    //Stupid monte-carlo integration of cos(x) from 0 to pi
    integrate_me(path, name.str(), node, n_nodes);
    
    //Blocking parameters
    int max = 1000;
    int min = n_nodes;
    int nblocks = 100;

    //start blocking. Process [i] reads a unique filename [name]
    Blocking* blocker = new Blocking(path, name.str(), max, min, nblocks, node, n_nodes);
    double finalErr = blocker->estimate_error();
    
    if (node == 0) cout << "Final stddev: " << finalErr << endl;
    
    MPI_Finalize();

}
