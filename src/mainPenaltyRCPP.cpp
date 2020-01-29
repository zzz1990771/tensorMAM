//[[Rcpp::depends(RcppEigen)]]
#include "mam.h"

//----------------------------------------------------------------**
//***--------------------Estimation without penalty---------------**
// [[Rcpp::export]]
List Estimation(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S, double threshold, int max_step)
{
	double  likhd0 = pow(10, 6), likhd1 = 0;
	MatrixXd Dnew;
	VectorXi convergence1;	
	
	int step = 0;
	while (step<max_step) {
		convergence1 = VectorXi::Constant(4, 1);
		step = step + 1;
		MatrixXd Snew = updateS(Y, Z, A, B, C);
		Dnew = C * Snew*kroneckerProduct(B.transpose(), A.transpose());
		likhd1 = (Y - Z * Dnew.transpose()).squaredNorm();
		if (likhd1<likhd0) {
			S = Snew;
			likhd0 = likhd1;
		}
		else convergence1[0]=0;
		MatrixXd Cnew = updateC(Y, Z, A, B, C, S);
		Dnew = Cnew * S*kroneckerProduct(B.transpose(), A.transpose());
		likhd1 = (Y - Z * Dnew.transpose()).squaredNorm();

		if (likhd1<likhd0) {
			C = Cnew;
			likhd0 = likhd1;
		}
		else convergence1[1]=0;

		MatrixXd Anew = updateA(Y, Z, A, B, C, S);
		Dnew = C * S*kroneckerProduct(B.transpose(), Anew.transpose());
		likhd1 = (Y - Z * Dnew.transpose()).squaredNorm();
		if (likhd1<likhd0) {
			A = Anew;
			likhd0 = likhd1;
		}
		else convergence1[2]=0;

		MatrixXd Bnew = updateB(Y, Z, A, B, C, S);
		Dnew = C * S*kroneckerProduct(Bnew.transpose(), A.transpose());
		likhd1 = (Y  - Z * Dnew.transpose()).squaredNorm();
		if (likhd1<likhd0) {
			B = Bnew;
			if ((likhd0 - likhd1) / likhd0<threshold) break;
			else  likhd0 = likhd1;		
		}
		else convergence1[3]=0;
		if(convergence1.sum()==0) break;
	}
	return List::create(Named("likhd") = likhd0, Named("Dnew") = Dnew,Named("S") = S);
}
//----------------------------------------------------------------**
//***--------------------setup tuning parameters------------------**
// [[Rcpp::export]]
VectorXd setuplambda(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S, int nlam, VectorXd setlam)
{
	int n = Y.rows(), q = Y.cols(), p = A.rows(), r1 = A.cols(), r2 = B.cols(), r3 = C.cols(), k = B.rows(), j, jj;

	double lam_max, lam_min, alpha;
	VectorXd lambda, lambda1, tmp, id,tmp1;
	MatrixXd S1, cbs, V, V_1, Y1 = Y, Gammaj, Gamma_sqrt, svdu, svdd, U;
	Y1.resize(n*q, 1);
	S1 = TransferModalUnfoldings(S, 3, 1, r1, r2, r3);
	cbs = kroneckerProduct(C, B)*(S1.transpose());
	id = SEQ(1, p*k, p);
	tmp = VectorXd::Constant(p, 0);
	for (j = 0; j < p; j++)
	{
		V = submatrix_col(Z, id.array() + j)*(cbs.block(0, 0, k, r1));
		for (jj = 1; jj < q; jj++) {
			V_1 = submatrix_col(Z, id.array() + j)*(cbs.block(jj*k, 0, k, r1));
			V = rbind_rcpp(V, V_1);
		}
		Gammaj = V.transpose()*V/ n;
		
		
		JacobiSVD<MatrixXd> svd(Gammaj, ComputeThinU | ComputeThinV);
		svdu = svd.matrixU();
		svdd = (svd.singularValues()).asDiagonal();
		Gamma_sqrt = svdu * ((1 / (svdd.diagonal().array().sqrt())).matrix().asDiagonal())*(svdu.transpose());	
		tmp1 = Y1.transpose()*V * Gamma_sqrt;
		tmp[j]=tmp1.array().abs().sum();
	}
	lam_max = setlam[0];
	lam_min = setlam[1];
	alpha = setlam[2];

	double max_tmp;
	max_tmp = (tmp.array()).maxCoeff()/sqrt(n*q*q);
	double max_lam;
	max_lam = lam_max * max_tmp / alpha;
	if (lam_min == 0) {
		lambda1.setLinSpaced(nlam - 1, log(max_lam), log(0.0001*max_lam));
		lambda.setLinSpaced(nlam, 0, 0);
		lambda.segment(0, nlam - 1) = lambda1.array().exp().matrix();
	}
	else {
		lambda1.setLinSpaced(nlam, log(max_lam), log(lam_min*max_lam));
		lambda = lambda1.array().exp();
	}
	return lambda;
}


//***-------------------------------------------------------------**
//***-------update the jth row of matrix A with penalty-----------**
MatrixXd updateA_penalty(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S, MatrixXd beta,
	VectorXd &activeA, double lambda1,int nlam,int max_iter, double alpha, double gamma, double penalty,int dfmax, double eps)
{
/*
	Input:
	Y is n*q matrix
	Z is n*(K*p) matrix
	A is p*r1 matrix
	B is K*r2 matrix
	C is q*r3 matrix
	S is r3*(r1*r2) matrix
	penalty= 1(lasso),2(mcp), and 3(scad)
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	eps is error to control convergence
	nlam is the number of preset tuning parameters
	max_iter is the maxixum number of iteration
*/  
	int r1 = A.cols();
	int r2 = B.cols();
	int r3 = C.cols();
	int p = A.rows();
	int k = B.rows();
	int n = Y.rows();
	int q = C.rows();
	int j,jj, active;

	MatrixXd S1,S_3;
	VectorXd aj,ajnew,zj;

	S1 = S.row(0).transpose();
	S1.resize(r1, r2);
	for (j = 1; j < r3; j++) {
		S_3 = S.row(j).transpose();
		S_3.resize(r1, r2);
		S1 = cbind_rcpp(S1, S_3);
	}
	VectorXd id;
	id = SEQ(1, p*k, p);
	MatrixXd Vnew;
	Vnew = MatrixXd::Constant(n*q, r1*p, 0);
	MatrixXd Gamma_sqrtn;
	Gamma_sqrtn = MatrixXd::Constant(r1, r1*p, 0);
	MatrixXd cbs;
	cbs = kroneckerProduct(C, B)*(S1.transpose());
	MatrixXd A1 = A, IDEN, V, Gammaj, Gamma_sqrt, V_1;
	IDEN = MatrixXd::Identity(q, q);
	MatrixXd D3,L;
	D3 = C * S*kroneckerProduct(B.transpose(), A.transpose());

	for (j = 0; j < p; j++) {
		V = submatrix_col(Z, id.array() + j)*(cbs.block(0, 0, k, r1));
		for (jj = 1; jj < q; jj++) {
			V_1 = submatrix_col(Z, id.array() + j)*(cbs.block(jj*k, 0, k, r1));
			V = rbind_rcpp(V, V_1);
		}
		Gammaj = ((V.transpose()*V).array() / n).matrix();
		L = Gammaj.llt().matrixL();
		Gamma_sqrtn.block(0, j*r1, r1, r1) = UpTriangularInv(L.transpose());
		Vnew.block(0, j*r1, n*q, r1) = QbyR(V, UpTriangularInv(L.transpose()),1);
		A1.row(j) = (QbyR(A.row(j).transpose(), L.transpose(), 0)).transpose();
	}		

	MatrixXd Anew = A1;
	MatrixXd r = Y - Z * (D3.transpose());
	r.resize(n*q, 1);
	VectorXd ajnorm_old, ajnorm;
	ajnorm_old = ajnorm = VectorXd::Constant(p, 0);
	double converged1;
	int step = 0;
	while (step<max_iter) 
	{
		step++;
		active = 0;
		for (j = 0; j < p; j++)
			if (ajnorm[j] != 0) active = active + 1;
		if (active>dfmax) {
			beta = MatrixXd::Constant(p*r1, nlam, -9);
			return A;
		}
		for (j = 0; j < p;j++) {
			aj = Anew.row(j).transpose();
			zj = Vnew.block(0, j*r1, n*q, r1).transpose()*r + n*aj;
			ajnew = updateAj(zj, n, r1, lambda1, alpha, gamma, penalty);
			r = r - Vnew.block(0, j*r1, n*q, r1)*(ajnew-aj);
			Anew.row(j) = ajnew.transpose();
			ajnorm[j] = ajnew.norm();
		}
		converged1 = 1;
		for (j = 0; j < p;j++) {
				if (ajnorm[j] != 0 && ajnorm_old[j] != 0) {
					if ((A1.row(j) - Anew.row(j)).norm() / ajnorm_old[j]>eps) {
						converged1 = 0; break;
					}
				}
				else if (ajnorm[j] == 0 && ajnorm_old[j] != 0) {
					converged1 = 0; break;
				}
				else if (ajnorm[j] != 0 && ajnorm_old[j] == 0) {
					converged1 = 0; break;
				}
			}
		if (converged1) break;
		A1 = Anew;
		ajnorm_old = ajnorm;
    }//end while
	for (j = 0; j<p; j++) {
		Anew.row(j) = (QbyR(Anew.row(j).transpose(),Gamma_sqrtn.block(0, j*r1, r1, r1),0)).transpose();
		if (ajnorm[j]) activeA[j] = 1;
	}
	return Anew;
}

//***-------------------------------------------------------------**
//***------Old main function: Estimation with penalizing functions in a whole column -----------------**
// [[Rcpp::export]]
List EstPenColumn(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S, VectorXd lambda,
	double alpha, double gamma, double penalty, int dfmax, double threshold, double eps, int max_step, int max_iter)
{
	/*
	Input:
	Y is n*q matrix
	Z is n*(K*p) matrix
	A is p*r1 matrix
	B is K*r2 matrix
	C is q*r3 matrix
	S is r3*(r1*r2) matrix
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	gamma is another tuning parameter for MCP and SCAD
	penalty= 1(lasso),2(mcp), and 3(scad)
	dfmax is the preset maximum digrees of freedom
	threshold is the error to control convergence for outer iteration
	eps is the error to control convergence for inner iteration
	max_step is the max step to control convergence for outer iteration
	max_iter is the max step to control convergence for inner iteration
	is_setlam is logical, 1 for set lambda by data; 0 for given lambda
	setlam is a vector saving preseting lam_max,lam_min,alpha,nlam
	nlam is the number of tuning parameters

	Output:
	Dnew is a estimator of D(3) = C%*%S%*%kronecker(t(B),t(A))
	*/
	int l,j, step, nlam, p = A.rows(), r1 = A.cols(), n = Y.rows(), K = B.rows();
	double  likhd0 = pow(10, 6), lambda1, likhd1 = 0;
	MatrixXd Dnew, Anew, Bnew, Snew, Cnew, Z1=Z, A1, Z2, A2;
	VectorXd activeA, convergence1;
	nlam = lambda.size();
	VectorXd likhd = VectorXd::Constant(nlam, 0);
	VectorXd df = VectorXd::Constant(nlam, 0);
	MatrixXd betapath;
	betapath = MatrixXd::Constant(p, nlam, 0);
	Anew = A1 = A;
	for (l = 0; l < nlam; l++) {
		lambda1 = n * lambda[l];
		step = 0;
		while (step<max_step) {
		  convergence1 = VectorXd::Constant(4, 1);
			step ++;
			Snew = updateS(Y, Z1, A1, B, C);
			Dnew = C * Snew * kroneckerProduct(B.transpose(), A1.transpose());
			likhd1 = (Y - Z1 * Dnew.transpose()).squaredNorm();
			if (likhd1<likhd0) {
				S = Snew;
				likhd0 = likhd1;
			}
			else convergence1[0]=0;
			Cnew = updateC(Y, Z1, A1, B, C, S);
			Dnew = Cnew * S * kroneckerProduct(B.transpose(), A1.transpose());
			likhd1 = (Y - Z1 * Dnew.transpose()).squaredNorm();
			if (likhd1<likhd0) {
				C = Cnew;
				likhd0 = likhd1;
			}
			else convergence1[1]=0;
			activeA = VectorXd::Constant(p, 0);
			Anew = updateA_penalty(Y, Z, A, B, C, S, betapath, activeA, lambda1, nlam, max_iter, alpha, gamma, penalty, dfmax, eps);
			if(activeA.sum()<r1){
				A1 = A;
				Z1 = Z;
				break;
			}
			else{
				Z2 = extractColsZ(Z,p,K,activeA);
				A2 = extractRows(Anew, activeA);
				Dnew = C * S * kroneckerProduct(B.transpose(), A2.transpose());
				likhd1 = (Y - Z2 * Dnew.transpose()).squaredNorm();
				if (likhd1<likhd0) {
					A = Anew;
					Z1 = Z2;
					A1 = A2;
					likhd0 = likhd1;
				}
				else convergence1[2]=0;
			}

			Bnew = updateB(Y, Z1, A1, B, C, S);
			Dnew = C * S * kroneckerProduct(Bnew.transpose(), A1.transpose());
			likhd1 = (Y - Z1 * Dnew.transpose()).squaredNorm();
			if (likhd1<likhd0) {
				B = Bnew;
				if ((likhd0 - likhd1) / likhd0<threshold) break;
				else  likhd0 = likhd1;	
			}
			else convergence1[3]=0;
		
			if(convergence1.sum()==0) break;
		} //end while		
		activeA = VectorXd::Constant(p, 0);
		for(j=0;j<p;j++) if(A.row(j).norm()) activeA[j] = 1;
		df[l] = activeA.sum();
		likhd[l] = likhd0;
		betapath.col(l) = activeA;
	}// end for
	Dnew = C * S * kroneckerProduct(B.transpose(), A.transpose());
	return List::create(Named("likhd") = likhd, Named("betapath") = betapath, Named("df") = df, Named("Dnew") = Dnew, Named("lambda")=lambda);
}

void updateA_penaltyl(VectorXd r, MatrixXd Ztilde, MatrixXd &Dl, MatrixXd beta, MatrixXd &activeA, int K, int p,
	int l, int r1, double lambda1, int nlam, int max_iter, double alpha, double gamma, double penalty, int dfmax, double eps)
{
	/*
	Input:
	r is n-vector, which is the lth column of residual Y - Z * (D3.transpose())
	Ztilde is n*Kp matrix, which is (\tilde{z}_1^{new},..., \tilde{z}_p^{new})
	Dl is K*p matrix, which is reshaped from lth row of D_{(3)}
	activeX is q*p matrix, which saves the active variables
	penalty= 1(lasso),2(mcp), and 3(scad)
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	eps is error to control convergence
	nlam is the number of preset tuning parameters
	max_iter is the maxixum number of iteration

	output:
	Dl is K*p matrix
*/
	int n = Ztilde.rows();
	int j, active;
	VectorXd djnorm_old, djnorm, dj, djnew, gj;
	MatrixXd Dlnew = Dl, Zj;
	djnorm_old = djnorm = VectorXd::Constant(p, 0);
	double converged1;
	int step = 0;
	
	while (step < max_iter)
	{
		step++;
		active = 0;
		for (j = 0; j < p; j++)
			if (djnorm[j] != 0) active = active + 1;
		if (active > dfmax) {
			beta = MatrixXd::Constant(p * r1, nlam, -9);
			break; 
		}
		for (j = 0; j < p; j++) {
			Zj = Ztilde.block(0, j * K, n, K);
			dj = Dlnew.col(j);
			gj = Zj.transpose() * r + n * dj;
			djnew = updateAj(gj, n, K, lambda1, alpha, gamma, penalty);
			r = r - Zj * (djnew - dj);
			Dlnew.col(j) = djnew;
			djnorm[j] = djnew.norm();
		}

		converged1 = 1;
		for (j = 0; j < p; j++) {
			if (djnorm[j] != 0 && djnorm_old[j] != 0) {
				if ((Dl.col(j) - Dlnew.col(j)).norm() / djnorm_old[j] > eps) {
					converged1 = 0; break;
				}
			}
			else if (djnorm[j] == 0 && djnorm_old[j] != 0) {
				converged1 = 0; break;
			}
			else if (djnorm[j] != 0 && djnorm_old[j] == 0) {
				converged1 = 0; break;
			}
		}
		if (converged1)  break;
		Dl = Dlnew;
		djnorm_old = djnorm;
	}//end while
	for (j = 0; j < p; j++) if (djnorm[j]) activeA(l, j) = 1;
}

MatrixXd updateA_penalty_single(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S, MatrixXd beta,
	MatrixXd &activeA, double lambda1, int nlam, int max_iter, double alpha, double gamma, double penalty, int dfmax, double eps)
{
	/*
		Input:
		Y is n*q matrix
		Z is n*(K*p) matrix
		A is p*r1 matrix
		B is K*r2 matrix
		C is q*r3 matrix
		S is r3*(r1*r2) matrix
		activeA is q*p matrix
		penalty= 1(lasso),2(mcp), and 3(scad)
		lambda is preset tuning parameters, a L-vector
		alpha is tuning parameter being in (0,1) to control elasnet
		eps is error to control convergence
		nlam is the number of preset tuning parameters
		max_iter is the maxixum number of iteration
	*/

	int i,j,r1 = A.cols(), r2 = B.cols(), r3 = C.cols(), p = A.rows(), K = B.rows(), n = Y.rows(), q = C.rows();
	MatrixXd S1, V, V_inv, Anew = A, R, Q;
	S1 = TransferModalUnfoldings(S, 3, 1, r1, r2, r3);
	V = kroneckerProduct(C, B) * (S1.transpose());
	HouseholderQR<MatrixXd> qr;
	qr.compute(V);
	Q = qr.householderQ();
	R = qr.matrixQR().triangularView<Upper>();
	V_inv = QbyR(Q.transpose(), UpTriangularInv(R), 0);


	VectorXd id, Rl;
	id = SEQ(1, p * K, p);
	MatrixXd Dl, Gamma_sqrtn, D3, Ztildej, D3tildej, Gammaj, r;
	Gamma_sqrtn = MatrixXd::Constant(K, K * p, 0);
	
	D3 = C * S * kroneckerProduct(B.transpose(), A.transpose());
	r = Y - Z * (D3.transpose());

	MatrixXd D3tilde_new_hat, ditilde, D3tilde_new, Ztilde_new, U;
	D3tilde_new_hat = D3tilde_new= MatrixXd::Constant(q, p * K, 0);
	Ztilde_new = MatrixXd::Constant(n, p * K, 0);	
	
	for (j = 0; j < p; j++) {
		Ztildej = submatrix_col(Z, id.array() + j);
		Gammaj = ((Ztildej.transpose() * Ztildej).array() / n).matrix();
		U = Gammaj.llt().matrixU();
		Gamma_sqrtn.block(0, j * K, K, K) = UpTriangularInv(U);
		
		D3tildej = submatrix_col(D3, id.array() + j);
		D3tilde_new.block(0, j * K, q, K) = QbyR(D3tildej.transpose(), U, 0).transpose();
		Ztilde_new.block(0, j * K, n, K) = QbyR(Ztildej, Gamma_sqrtn.block(0, j * K, K, K), 1);
	}
	for (i = 0; i < q; i++) {
		Rl = r.col(i);
		Dl = D3tilde_new.row(i).transpose();
		Dl.resize(K, p);
		updateA_penaltyl(Rl, Ztilde_new, Dl, beta, activeA, K, p, i, r1, lambda1, nlam, max_iter, alpha, gamma, penalty, dfmax, eps);
		Dl.resize(K * p, 1);
		D3tilde_new_hat.row(i) = Dl.transpose();	
	}	
	for (j = 0; j < p; j++) {
		ditilde = QbyR(D3tilde_new_hat.block(0, j * K, q, K).transpose(), Gamma_sqrtn.block(0, j * K, K, K), 0);
		ditilde.resize(K * q, 1);
		Anew.row(j) = (V_inv * ditilde).transpose();
	}
	return(Anew);
}

//***-------------------------------------------------------------**
//***------main function: Estimation with penalizing single function -----------------**
// [[Rcpp::export]]
List EstPenSingle(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S, VectorXd lambda,
	double alpha, double gamma, double penalty, int dfmax, double threshold, double eps, int max_step, int max_iter)
{
	/*
	Input:
	Y is n*q matrix
	Z is n*(K*p) matrix
	A is p*r1 matrix
	B is K*r2 matrix
	C is q*r3 matrix
	S is r3*(r1*r2) matrix
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	gamma is another tuning parameter for MCP and SCAD
	penalty= 1(lasso),2(mcp), and 3(scad)
	dfmax is the preset maximum digrees of freedom
	threshold is the error to control convergence for outer iteration
	eps is the error to control convergence for inner iteration
	max_step is the max step to control convergence for outer iteration
	max_iter is the max step to control convergence for inner iteration
	is_setlam is logical, 1 for set lambda by data; 0 for given lambda
	setlam is a vector saving preseting lam_max,lam_min,alpha,nlam
	nlam is the number of tuning parameters

	Output:
	Dnew is a estimator of D(3) = C%*%S%*%kronecker(t(B),t(A))
	*/
	int l, j, step, nlam, p = A.rows(), r1 = A.cols(), n = Y.rows(), q = C.rows();
	double  likhd0 = pow(10, 6), lambda1, likhd1 = 0;
	MatrixXd Dnew, Anew, Bnew, Snew, Cnew;	
	
	VectorXd convergence1, activeX;
	MatrixXd activeA;
	nlam = lambda.size();
	VectorXd df, likhd;
	likhd = df = VectorXd::Constant(nlam, 0);
	MatrixXd betapath, activeXpath;
	betapath = MatrixXd::Constant(q*p, nlam, 0);
	activeXpath = MatrixXd::Constant(p, nlam, 0);
	Anew = A;
	for (l = 0; l < nlam; l++) {
		lambda1 = n * lambda[l];
		step = 0;
		while (step<max_step) {
		  convergence1 = VectorXd::Constant(4, 1);
			step ++;
			Snew = updateS(Y, Z, A, B, C);
			Dnew = C * Snew * kroneckerProduct(B.transpose(), A.transpose());
			likhd1 = (Y - Z * Dnew.transpose()).squaredNorm();
			if (likhd1<likhd0) {
				S = Snew;
				likhd0 = likhd1;
			}
			else convergence1[0]=0;

			Cnew = updateC(Y, Z, A, B, C, S);
			Dnew = Cnew * S * kroneckerProduct(B.transpose(), A.transpose());
			likhd1 = (Y - Z * Dnew.transpose()).squaredNorm();
			if (likhd1<likhd0) {
				C = Cnew;
				likhd0 = likhd1;
			}
			else convergence1[1]=0;
			activeA = MatrixXd::Constant(q, p, 0);
			Anew = updateA_penalty_single(Y, Z, A, B, C, S, betapath, activeA, lambda1, nlam, max_iter, alpha, gamma, penalty, dfmax, eps);
			activeX = VectorXd::Constant(p, 0);
			for(j=0;j<p;j++) if(Anew.row(j).norm()) activeX[j]= 1;
			if(activeX.sum()<r1)	break;
			else{
				Dnew = C * S * kroneckerProduct(B.transpose(), Anew.transpose());
				likhd1 = (Y - Z * Dnew.transpose()).squaredNorm();
				if (likhd1<likhd0) {
					A = Anew;
					likhd0 = likhd1;
				}
				else convergence1[2]=0;
			}
			Bnew = updateB(Y, Z, A, B, C, S);
			Dnew = C * S * kroneckerProduct(Bnew.transpose(), A.transpose());
			likhd1 = (Y - Z * Dnew.transpose()).squaredNorm();
			if (likhd1<likhd0) {
				B = Bnew;			
				if ((likhd0 - likhd1) / likhd0<threshold) break;
				else  likhd0 = likhd1;	
			}
			else convergence1[3]=0;		
			if(convergence1.sum()==0) break;
		} //end while
		activeX = VectorXd::Constant(p, 0);
		for(j=0;j<p;j++) if(A.row(j).norm()) activeX[j]= 1;		
		df[l] = activeX.sum();
		likhd[l] = likhd0;
		activeA.resize(q*p,1);
		betapath.col(l) = activeA;
		activeXpath.col(l) = activeX;
	}// end for		
	Dnew = C * S * kroneckerProduct(B.transpose(), A.transpose());	
	return List::create(Named("likhd") = likhd, Named("betapath") = betapath, Named("df") = df, Named("Dnew") = Dnew, 
	                    Named("lambda")=lambda, Named("activeXpath") = activeXpath);
}

//***--------------------------------------------------------------**
//***----------- Estimation with penalizing functions in a whole column by CV---------------------**
// [[Rcpp::export]]
List EstPenColumnCV(MatrixXd Y, MatrixXd Z, MatrixXd Ytest, MatrixXd Ztest, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S, VectorXd lambda,
	double alpha, double gamma, double penalty, int dfmax, double threshold, double eps, int max_step, int max_iter)
{
	/*
	Input:
	Y is n*q matrix
	Z is n*(K*p) matrix
	A is p*r1 matrix
	B is K*r2 matrix
	C is q*r3 matrix
	S is r3*(r1*r2) matrix
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	gamma is another tuning parameter for MCP and SCAD
	penalty= 1(lasso),2(mcp), and 3(scad)
	dfmax is the preset maximum digrees of freedom
	threshold is the error to control convergence for outer iteration
	eps is the error to control convergence for inner iteration
	max_step is the max step to control convergence for outer iteration
	max_iter is the max step to control convergence for inner iteration
	is_setlam is logical, 1 for set lambda by data; 0 for given lambda
	setlam is a vector saving preseting lam_max,lam_min,alpha,nlam
	nlam is the number of tuning parameters

	Output:
	Dnew is a estimator of D(3) = C%*%S%*%kronecker(t(B),t(A))
	*/
	int l,j, step, nlam, p = A.rows(), r1 = A.cols(), n = Y.rows(), K = B.rows();
	double  likhd0 = pow(10, 6), lambda1, likhd1 = 0;
	MatrixXd Dnew, Anew, Bnew, Snew, Cnew, Z1=Z, A1, Z2, A2;	
	VectorXd activeA, convergence1;
	nlam = lambda.size();
	VectorXd likhd = VectorXd::Constant(nlam, 0);
	VectorXd df = VectorXd::Constant(nlam, 0);
	MatrixXd betapath = MatrixXd::Constant(p*r1, nlam, 0);
	Anew = A1 = A;
	for (l = 0; l < nlam; l++) {
		lambda1 = n * lambda[l];
		step = 0;
		while (step<max_step) {
		  convergence1 = VectorXd::Constant(4, 1);
			step ++;
			Snew = updateS(Y, Z1, A1, B, C);
			Dnew = C * Snew * kroneckerProduct(B.transpose(), A1.transpose());
			likhd1 = (Y - Z1 * Dnew.transpose()).squaredNorm();
			if (likhd1<likhd0) {
				S = Snew;
				likhd0 = likhd1;
			}
			else convergence1[0]=0;
			Cnew = updateC(Y, Z1, A1, B, C, S);
			Dnew = Cnew * S * kroneckerProduct(B.transpose(), A1.transpose());
			likhd1 = (Y - Z1 * Dnew.transpose()).squaredNorm();
			if (likhd1<likhd0) {
				C = Cnew;
				likhd0 = likhd1;
			}
			else convergence1[1]=0;
			activeA = VectorXd::Constant(p, 0);
			Anew = updateA_penalty(Y, Z, A, B, C, S, betapath, activeA, lambda1, nlam, max_iter, alpha, gamma, penalty, dfmax, eps);
			if(activeA.sum()<r1){
				A1 = A;
				Z1 = Z;
				break;
			}
			else{
				Z2 = extractColsZ(Z,p,K,activeA);
				A2 = extractRows(Anew, activeA);
				Dnew = C * S * kroneckerProduct(B.transpose(), A2.transpose());
				likhd1 = (Y - Z2 * Dnew.transpose()).squaredNorm();
				if (likhd1<likhd0) {
					A = Anew;
					Z1 = Z2;
					A1 = A2;
					likhd0 = likhd1;
				}
				else convergence1[2]=0;
			}
			Bnew = updateB(Y, Z1, A1, B, C, S);
			Dnew = C * S * kroneckerProduct(Bnew.transpose(), A1.transpose());
			likhd1 = (Y - Z1 * Dnew.transpose()).squaredNorm();
			if (likhd1<likhd0) {
				B = Bnew;				
				if ((likhd0 - likhd1) / likhd0<threshold) break;
				else  likhd0 = likhd1;				
			}
			else convergence1[3]=0;		
			if(convergence1.sum()==0) break;
		} //end while
		activeA = VectorXd::Constant(p, 0);
		for(j=0;j<p;j++) if(A.row(j).norm()) activeA[j] = 1;
		df[l] = activeA.sum();
		Z2 = extractColsZ(Ztest,p,K,activeA);
		A2 = extractRows(A, activeA);
		Dnew = C * S * kroneckerProduct(B.transpose(), A2.transpose());
		likhd[l] = (Y - Z2 * Dnew.transpose()).squaredNorm();
		betapath.col(l) = activeA;
	}// end for
	return List::create(Named("likhd") = likhd, Named("df") = df, Named("betapath")=betapath);
}

//***--------------------------------------------------------------**
//***----------- Estimation with penalizing single function by CV---------------------**
// [[Rcpp::export]]
List EstPenSingleCV(MatrixXd Y, MatrixXd Z, MatrixXd Ytest, MatrixXd Ztest, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S, VectorXd lambda,
	double alpha, double gamma, double penalty, int dfmax, double threshold, double eps, int max_step, int max_iter)
{
	/*
	Input:
	Y is n*q matrix
	Z is n*(K*p) matrix
	A is p*r1 matrix
	B is K*r2 matrix
	C is q*r3 matrix
	S is r3*(r1*r2) matrix
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	gamma is another tuning parameter for MCP and SCAD
	penalty= 1(lasso),2(mcp), and 3(scad)
	dfmax is the preset maximum digrees of freedom
	threshold is the error to control convergence for outer iteration
	eps is the error to control convergence for inner iteration
	max_step is the max step to control convergence for outer iteration
	max_iter is the max step to control convergence for inner iteration
	is_setlam is logical, 1 for set lambda by data; 0 for given lambda
	setlam is a vector saving preseting lam_max,lam_min,alpha,nlam
	nlam is the number of tuning parameters

	Output:
	Dnew is a estimator of D(3) = C%*%S%*%kronecker(t(B),t(A))
	*/
	int l,j, step, nlam, p = A.rows(), r1 = A.cols(), n = Y.rows(), q = C.rows();
	double  likhd0 = pow(10, 6), lambda1, likhd1 = 0;
	MatrixXd Dnew, Anew, Bnew, Snew, Cnew, activeA;
	VectorXd convergence1, activeX, likhd, df;
	nlam = lambda.size();
	likhd = df = VectorXd::Constant(nlam, 0);	
	MatrixXd betapath, activeXpath;
	betapath = MatrixXd::Constant(q*p, nlam, 0);
	activeXpath = MatrixXd::Constant(p, nlam, 0);	
	Anew = A;
	for (l = 0; l < nlam; l++) {
		lambda1 = n * lambda[l];
		step = 0;
		while (step<max_step) {
		  convergence1 = VectorXd::Constant(4, 1);
			step ++;
			Snew = updateS(Y, Z, A, B, C);
			Dnew = C * Snew * kroneckerProduct(B.transpose(), A.transpose());
			likhd1 = (Y - Z * Dnew.transpose()).squaredNorm();
			if (likhd1<likhd0) {
				S = Snew;
				likhd0 = likhd1;
			}
			else convergence1[0]=0;
			Cnew = updateC(Y, Z, A, B, C, S);
			Dnew = Cnew * S * kroneckerProduct(B.transpose(), A.transpose());
			likhd1 = (Y - Z * Dnew.transpose()).squaredNorm();
			if (likhd1<likhd0) {
				C = Cnew;
				likhd0 = likhd1;
			}
			else convergence1[1]=0;
			activeA = MatrixXd::Constant(q, p, 0);
			Anew = updateA_penalty_single(Y, Z, A, B, C, S, betapath, activeA, lambda1, nlam, max_iter, alpha, gamma, penalty, dfmax, eps);
			activeX = VectorXd::Constant(p, 0);
			for(j=0;j<p;j++) if(Anew.row(j).norm()) activeX[j]= 1;			
			if(activeX.sum()<r1)  break;
			else{
				Dnew = C * S * kroneckerProduct(B.transpose(), Anew.transpose());
				likhd1 = (Y - Z * Dnew.transpose()).squaredNorm();
				if (likhd1<likhd0) {
					A = Anew;
					likhd0 = likhd1;
				}
				else convergence1[2]=0;
			}
			Bnew = updateB(Y, Z, A, B, C, S);
			Dnew = C * S * kroneckerProduct(Bnew.transpose(), A.transpose());
			likhd1 = (Y - Z * Dnew.transpose()).squaredNorm();		
			if (likhd1<likhd0) {
				B = Bnew;				
				if ((likhd0 - likhd1) / likhd0<threshold) break;
				else  likhd0 = likhd1;				
			}
			else convergence1[3]=0;		
			if(convergence1.sum()==0) break;
		} //end while		
		activeX = VectorXd::Constant(p, 0);
		for(j=0;j<p;j++) if(A.row(j).norm()) activeX[j]= 1;		
		df[l] = activeX.sum();		
		Dnew = C * S * kroneckerProduct(B.transpose(), A.transpose());		
		likhd[l] = (Y - Ztest * Dnew.transpose()).squaredNorm();				
		activeA.resize(q*p,1);
		betapath.col(l) = activeA;
		activeXpath.col(l) = activeX;		
	}// end for		
	return List::create(Named("likhd") = likhd, Named("df") = df, Named("betapath")=betapath, Named("activeXpath") = activeXpath);
}

//----------------------------------------------------------------**
//***--------------------EstimationD3 directly--------------------**
// [[Rcpp::export]]
List EstimationD3(MatrixXd Y, MatrixXd Z)
{
  int j,q = Y.cols(),kp = Z.cols();
  MatrixXd Dnew = MatrixXd::Constant(q, kp, 0);
  HouseholderQR<MatrixXd> qr;
  qr.compute(Z);
  MatrixXd R = qr.matrixQR().triangularView<Upper>();
  MatrixXd Q = qr.householderQ();
  MatrixXd RQ = UpTriangularInv(R) * Q.transpose();
  
  MatrixXd temp, IDEN = MatrixXd::Identity(kp, kp);
  if (pow(condition_numberQRSym(R),2) > 1e10){
    temp = tRbyR(R) + (IDEN.array()*1e-4).matrix();
    for (j = 0; j < q; j++) Dnew.row(j) = (temp.colPivHouseholderQr().solve(Z.transpose()*Y.col(j))).transpose();
  }
  else
    for (j = 0; j < q; j++) Dnew.row(j) = (RQ * Y.col(j)).transpose();
  double likhd = (Y - Z * Dnew.transpose()).squaredNorm();
  return List::create(Named("likhd") = likhd, Named("Dnew") = Dnew);
}	
