#include <iostream>
#include <set>

//bio++
#include <Bpp/Numeric/Random.all>
//#include <NumCalc/random>
#include <Bpp/Numeric/RandomTools.h>
//#include <NumCalc/RandomTools.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/io.all>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTools.h>

// BOOST ublas
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/blas.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/mpi.hpp>


// BOOST ublas blas (& lapack, atlas etc.) bindings, nonstandard
// cf. http://svn.boost.org/svn/boost/sandbox/numeric_bindings/boost/numeric/bindings/
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/blas/blas.hpp>


typedef long double scalar_type;

typedef boost::numeric::ublas::vector<short> int_vector_type;

typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::vector<long double> long_vector_type;
typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> array_type;
typedef bpp::TreeTemplate<bpp::Node> tree_type;
typedef bpp::Node node_type;
typedef std::vector < std::vector <std::string> > scatter_string;
typedef std::vector < std::vector <scalar_type> > scatter_scalar;


//#mpi time_order
std::pair<Species_tree *,std::string > init_LL_mpi(std::string tree_file,std::string forest_file, int outgroup
		 ,const boost::mpi::communicator  world);
std::pair<Species_tree *,std::string > sim_init_LL_mpi(std::string tree_file,std::string forest_file, int outgroup
		 ,const boost::mpi::communicator  world);

scalar_type LL_mpi(Species_tree * infer_tree, std::string Sstring ,std::string mode,bool outgroup=true,scalar_type delta=0.01,scalar_type tau=0.01, scalar_type lambda=0.01);

scalar_type sum_counts(Species_tree * infer_tree,  const boost::mpi::communicator  world);
std::string strip_out(std::string Sstring);

void sum_ttf(Species_tree * infer_tree,  const boost::mpi::communicator  world);
std::string strip_out(std::string Sstring);

std::pair<scalar_type,std::string> SPR_step(scalar_type pre_ll,Species_tree * infer_tree, std::string Sstring,std::string mode="clock",bool greedy=true,bool root_step=false, bool skip_step=false);
