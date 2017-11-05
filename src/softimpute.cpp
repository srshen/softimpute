#include<RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


//[[Rcpp::export]]
arma::vec cutoff(arma::vec& s_cu, const double& lam_cu){
    arma::vec s_cut = zeros(size(s_cu));
    for(int i=0; i<s_cu.n_elem; i++){
        s_cut(i) = s_cu(i) - lam_cu;
        if(s_cut(i) < 0) s_cut(i) = 0;
    }
    return(s_cut);
}


//[[Rcpp::export]]
arma::mat S_lambda(const arma::mat& X_sl, const double& lam_sl){
    arma::mat U, V;
    arma::vec s;
    arma::mat Sigma = zeros(size(X_sl));
    svd(U,s,V,X_sl);
    arma::vec s_cutoff = cutoff(s,lam_sl);
    for(int i=0; i<s.n_elem;i++){
        Sigma(i,i) = s_cutoff(i);
    }
    arma::mat S_lam = U * Sigma * V.t();
    return S_lam; //only used for square matrix
}

//[[Rcpp::export]]
arma::mat P_Ome(const arma::mat X_p, const arma::mat P_O){
    arma::mat X_po = X_p;
    for(int i=0;i<X_p.n_rows;i++){
        for(int j=0;j<X_p.n_cols;j++){
            if(P_O(i,j)==0) X_po(i,j) = 0;
        }
    }
    return X_po;
}

//[[Rcpp::export]]
Rcpp::List softimpute(const arma::mat& X_miss, const arma::mat& P_si, const double& lambda){
    
    
    int nrow = X_miss.n_rows;
    int ncol = X_miss.n_cols;
    int maxiter = 200;
    double epsilon = 1e-5;
    
    arma::mat Z_old = zeros(nrow, ncol);
    arma::mat One   = ones(nrow, ncol);
    
    arma::mat temp  = P_Ome(X_miss, P_si) + P_Ome(Z_old, One - P_si);
    arma::mat Z_new = S_lambda(temp, lambda);
    
    int count = 0;
    double tune_par = 1;
    while(tune_par >= epsilon && count<maxiter){
        count ++;
        Z_old = Z_new + 0.5*(Z_new - Z_old); //accelerated.
        temp  = P_Ome(X_miss, P_si) + P_Ome(Z_old, One - P_si);
        Z_new = S_lambda(temp, lambda);
        for(int i=0;i<nrow;i++){
            for(int j=0;j<ncol;j++){
                if(Z_new(i,j)<0)  Z_new(i,j)=0;
                if(Z_new(i,j)>1)  Z_new(i,j)=1;
            }
        }
        tune_par = trace((Z_new - Z_old).t() * (Z_new - Z_old)) / trace(Z_old.t() * Z_old);
    }
    
    
    return List::create(Named("num_iteration") = count,
                        Named("Z")             = Z_new);
}


