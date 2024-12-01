#include <bits/stdc++.h>
using namespace std;
#define endl "\n"
#define int long long
#define vi vector<int>
#define vb vector<bool>
#define vd vector<double>
#define vc vector<char>
#define db double
#define si set<int>
#define pb push_back
#define vvi vector<vi>
#define vvb vector<vb>
#define vvd vector<vd>
#define vvc vector<vc>
#define vs vector<string>

const int N = 10;
const db h = 1/db(N);
const db q = 7/(4*M_PI);
const db beta = 0;
const db mu = 1;
const db T0 = 1;
const db C = 0.3;
const db s0 = 0.5;
const db rho0 =1;
const db u1 = 0.5;
const db xip = 0.5;
const db xis0 = 0.5;
const db xin0 = 1-xis0;
const db nu = 0.5;

vd prodmatvec(vvd& mat, vd& v) {
    // Check if multiplication is possible
    if (vec.size() != mat.size()) {
        cerr << "Error: Vector size and matrix row size must match!" << endl;
        return 1;
    }
    // Resulting vector size will be equal to the number of columns in the matrix
    vector<double> result(mat[0].size(), 0.0);

    // Perform vector-matrix multiplication
    for (size_t i = 0; i < mat[0].size(); ++i) {       // Loop through columns of the matrix
        for (size_t j = 0; j < vec.size(); ++j) {      // Loop through rows of the matrix/vector
            result[i] += vec[j] * mat[j][i];
        }
    }

    return result;
}

vd vectorsum(const vd& v1, const vd& v2) {
    if (v1.size() != v2.size()) {
        cout << "Vectors must be of the same size!" << endl;
        stopProgram();
    }
    int k = v1.size()

    vd v(k);
    for (int i = 0; i< k; ++i) {
        v[i] = v1[i] + v2[i];
    }
    
    return v;
}

vd vectordif(const vd& v1, const vd& v2) {
    if (v1.size() != v2.size()) {
        cout << "Vectors must be of the same size!" << endl;
        stopProgram();
    }
    int k = v1.size()

    vd v(k);
    for (int i = 0; i< k; ++i) {
        v[i] = v1[i] - v2[i];
    }
    
    return v;
}

vd vectorscale(const vd& v, db l) {
    int k = v.size();
    vd b(k);
    for (int i = 0; i < k; ++i) {
        b[i] = v[i] * l;
    }
    return b;
}

vvd matrixscale(const vvd& M, db l) {
    vvd result(M.size(), M[0].size());
    for (int i = 0; i < M.size(); i++) {
        for (int j = 0; j< M[0].size(); ++j) {
            result[i][j] = M[i][j]*l;
        }
    }
    return result;
}

vvd matrixsum(const vvd& A, const vvd& B) {
    vvd result(A.size(), A[0].size());
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j< A[0].size(); ++j) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result;
}

double dot(vd& v1, vd& v2) {
    return v1[0]*v2[0] + v1[1]*v2[1];
}

double norm(const vd& v) {
    db s = 0;
    int k = v.size();
    for (int i = 0; i <k; ++i) {
        s += v[i]*v[i];
    }
    return sqrt(s);
}

db w(db r) {
    if (r > 2*h) return 0;
    else return q/pow(h, 2)*pow(1-r/(2*h), 4)*(1+n*r/h);
}

db wp(db r) {
    if (r < 2*h) return 0;
    else return -5*r*q/(pow(h, 4)) * pow(1-r/(2*h), 3); //This formula may be wrong
}

vd rhodot(vvd& r, vvd& v, vd& m) {
    vd rhodot(N);
    for (int a = 0; a< N; ++a) {
        db Sa = 0;
        for (int b = 0; b < N; ++b) {
            if (b != a) {
                vd rab = vectordif(r[a], r[b]);
                vd vab = vectordif(v[a], v[b]);
                db wpab = wp(norm(rab));
                Sa += m[b]*wpab / norm(rab) * dot(rab, vab)
            }
        }
        rhodot[a] = Sa;
    }

    return rhodot;
}

db T(db s) {
    return (1 + (s-s0)/C)*T0;
}

db p(db rho) {
    return u1*u1*(rho-rho0);
}

db xis(db s) {
    return xis0 - xip*(s-s0);
}

db xin(db s) {
    return 1 - xis(s);
}

vd vn(vd v, vd vs, db s) {
    return vectorscale(vectordif(v, vectorscale(vs, xis(s))), 1/xin(s));
} 

vd vns(vd v, vd vs, db s) {
    return vectorscale(vectordif(v, vs), 1/xin(s));
}

vd theta(vd& m, vd& rho, vvd& r, vd& s, vvd& v, vvd& vs) {
    vd th(N);
    for (int a = 0; a < N; ++a) {
        db t1 = -beta;
        db t2 = -8*mu;
        db s1 = 0;
        db s2 = 0;

        for (int b = 0; b < N; ++b) {
            if (b != a) {
                vd rab = vectordif(r[a], r[b]);
                db x = norm(rab);
                db wpab = wp(x);
                s1 += m[b]/rho[b]*wpab/x*pow(T(s[a])-T(s[b]), 2);
                s2 += m[b]/rho[b]*(pow(dot(vectordif(vn(v[a], vs[a], s[a]), vn(v[b], vs[b], s[b])), 
                                            vectordif(r[a], r[b])), 2)) / (x*x + nu*nu) * wpab(x) / x;
            }
        }

        th[a] =  t1*s1 + t2*s2;
    }
    return th;
}

vd j(db rho,db& s, db& xis, vd& vs, vd& v) {
    return vectorscale(vns(v, vs, s), rho*s*xis);
}

vd sdot(vd& m, vd& rho, vvd& r, vd& s, vvd& v, vvd& vs) {
    vd Sdot(N);

    vd th = theta(m, rho, r, s, v, vs);
    vd Temp(N);
    vd Xis(N);
    vvd J(N);
    for (int a = 0; a < N; ++a) {
        Temp[a] = T(s[a]);
        Xis[a] = xis(s[a]);
        J[a] = j(rho[a], s[a], Xis[a], vs[a], v[a]);
    }
    
    for (int a = 0; a < N; ++a) {
        Sdot[a] = th[a]/ T(s[a]) / rho[a];
        db xisa = xis(s[a]);
        vd vnsa = vns(v[a], vs[a], s[a]);
        vd ja = j(rho[a], s[a], xisa, vnsa)
        for (int b = 0; b < N; ++b) {
            if (b != a) {
                vd rab = vectordif(r[a], r[b]);
                db x = norm(rab);
                db wpab = wp(x);
                
                Sdot[a] -= m[b]*wpab/x*(1/(rho[a]*rho[a]) * dot(J[a], rab)
                            + 1/(rho[b]*rho[b]) * dot(J[b], rab) 
                            - 2*beta*(Temp[a]-Temp[b])/rho[a]/rho[b]);
            }
        }

    }
    return Sdot;
}

vvd pi(db p, db rho, db xin, db xis, vd vns) {
    vvd PI(2, vd(2));
    PI[0][0] = p;
    PI[1][1] = p;

    db k = rho*xis*xin;
    for (int i = 0; i< 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            PI[i][j] += vns[i]*vns[j]*k;
        }
    }

    return PI;
}

vvd vdot(vd& m, vd& rho, vvd& r, vd& s, vvd& v, vvd& vs) {
    vd Vdot(N);
    vd Xin(N);
    vd Xis(N);
    vvd Vns(N);
    vector <vvd> Pi(N);
    vvd Vn(N);
    vd pres(N);
    for (int a = 0; a < N; ++a) {
        Xin[a] = xin(s[a]);
        Xis[a] = 1-Xin[a];
        Vns[a] = vns(v[a], vs[a], s[a]);
        pres[a] = p(rho[a]);
        Pi[a] = pi(pres[a], rho[a], Xin[a], Xis[a], Vns[a]);
        Vn[a] = vn(v[a], vs[a], s[a]);
    }

    for (int a = 0; a < N; ++a) {
        vd sum = {0, 0};
        for (int b = 0; b < N; ++b) {
            if (b != a) {
                vd rab = vectordif(r[a], r[b]);
                db x = norm(rab);
                db wpab = wp(x);

                vectordif(sum, vectorscale(prodmatvec(matrixsum(matrixscale(Pi[a], 1/(rho[a]*rho[a])), 
                                                            matrixscale(Pi[b], 1/(rho[b]*rho[b]))), rab) ,
                                             m[b]*wpab/x));
                                             
                vectorsum(sum, )
            }
        }
    }



}



void simulation() {
    int N = 10;
    int l = 1;

    vvd r(N, vd(2));
    vvd v(N, vd(2));
    vvd vs(N, vd(2));
    vd rho(N);
    vd s(N);
    vd T(N);
    vd m(N);



    


}

signed main(){
    int t;
    cin >> t;
    for (int i = 0; i < t; ++i) {
        
        solution();
    }
}
