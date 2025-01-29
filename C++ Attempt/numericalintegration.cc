#include <bits/stdc++.h>
using namespace std;
#define endl "\n"
#define vi vector<int>
#define vb vector<bool>
#define vd vector<double>
#define db double
#define pb push_back
#define vvi vector<vi>
#define vvd vector<vd>

const int N = 50;
const db q = 7/(4*M_PI);
const db beta = 0;
const db mu = 0;
const db T0 = 1.9;
const db C = 3902;
const db s0 = 0.5;
const db rho0 =145.5;
const db u1 = 40;
const db xip = 18.83;
const db xin0 = 0.4195;
const db xis0 = 1-xin0;
const db nu = 0.05;     //Non specified
const db dt = 0.00001;
const int n = 201;
const db L = 1;
const db h = 3*L/N;

vd prodmatvec(const vvd& mat, const vd& vec) {
    // Check if multiplication is possible
    if (vec.size() != mat.size()) {
        cerr << "Error: Vector size and matrix row size must match!" << endl;
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
    }
    int k = v1.size();

    vd v(k);
    for (int i = 0; i< k; ++i) {
        v[i] = v1[i] + v2[i];
    }
    
    return v;
}

vd vectordif(const vd& v1, const vd& v2) {
    if (v1.size() != v2.size()) {
        cout << "Vectors must be of the same size!" << endl;
    }
    int k = v1.size();

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
    int m = M.size();
    int n = M[0].size();
    vvd result(m, vd(n));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = M[i][j]*l;
        }
    }
    return result;
}

vvd matrixsum(const vvd& A, const vvd& B) {
    int m = A.size();
    int n = A[0].size();
    vvd result(m, vd(n));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j< n; ++j) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result;
}

db dot(const vd& v1, const vd& v2) {
    int k = v1.size();
    int m = v2.size();
    if (k != m) {
        cout << "Vectors must be of the same size!" << endl;
    }
    db sum = 0;
    for (int i = 0; i < k; ++i) {
        sum += v1[i]*v2[i];
    }
    return sum;
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
    else return q/pow(h, 2)*pow(1-r/(2*h), 4)*(1+2*r/h);
}

db wp(db r) {
    if (r < 2*h) return 0;
    else return -5*r*q/(pow(h, 4)) * pow(1-r/(2*h), 3); //This formula may be wrong
}

// The paper uses a closed formula
vd rhodot(vvd& r, vvd& v, vd& m) {
    vd rhodot(N);
    for (int a = 0; a< N; ++a) {
        db Sa = 0;
        for (int b = 0; b < N; ++b) {
            if (b != a) {
                vd rab = vectordif(r[a], r[b]);
                vd vab = vectordif(v[a], v[b]);
                db wpab = wp(norm(rab));
                Sa += m[b]*wpab / norm(rab) * dot(rab, vab);
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
                                            vectordif(r[a], r[b])), 2)) / (x*x + nu*nu) * wpab / x;
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
    vvd Vdot(N);
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
        vd sum(2, 0.0);
        for (int b = 0; b < N; ++b) {
            if (b != a) {
                vd rab = vectordif(r[a], r[b]);
                db x = norm(rab);
                db wpab = wp(x);
                vd vnab = vectordif(Vn[a], Vn[b]);

                sum = vectordif(sum, vectorscale(prodmatvec(matrixsum(matrixscale(Pi[a], 1/(rho[a]*rho[a])), 
                                                            matrixscale(Pi[b], 1/(rho[b]*rho[b]))), rab) ,
                                             m[b]*wpab/x));
                                             
                sum = vectorsum(sum, vectorscale(rab, m[b]*wpab/x*2*4*mu/rho[a]/rho[b]/(x*x + nu*nu)*dot(vnab, rab)));
            }
        }
        Vdot[a] = sum;
    }
    return Vdot;
}

vvd vsdot(vd& m, vd& rho, vvd& r, vd& s, vvd& v, vvd& vs) {
    vvd Vsdot(N);
    vd Xin(N);
    vd Xis(N);
    vvd Vns(N);
    vvd Vn(N);
    vd pres(N);
    vd Temp(N);
    for (int a = 0; a < N; ++a) {
        Xin[a] = xin(s[a]);
        Xis[a] = 1-Xin[a];
        Vns[a] = vns(v[a], vs[a], s[a]);
        pres[a] = p(rho[a]);
        Vn[a] = vn(v[a], vs[a], s[a]);
        Temp[a] = T(s[a]);
    }

    for (int a = 0; a < N; ++a) {
        vd sum(2, 0);
        for (int b = 0; b < N; ++b) {
            vd rab = vectordif(r[a], r[b]);
            db x = norm(rab);
            db wpab = wp(x);
            vd vnab = vectordif(Vn[a], Vn[b]);            
            if (b != a) {
                sum = vectordif(sum, vectorscale(rab, m[b]*wpab/x*Xin[a]/rho[a]*dot(vnab, Vns[a])));
                sum = vectordif(sum, vectorscale(rab, m[b]*wpab/x*(pres[a]/(rho[a]*rho[a]) + pres[b]/(rho[b]*rho[b]) + s[a]/rho[a]*(Temp[a] - Temp[b]))));
            }
        }
        Vsdot[a] = sum;
    }
    return Vsdot;
}

vd rho(vd& m, vvd& r, vd& C) {
    vd RHO = C;
    for (int a = 0;  a < N; ++a) {
        for (int b = 0; b < N; ++b) {
            if (b != a) {
                vd rab = vectordif(r[a], r[b]);
                db x = norm(rab);
                db wab = w(x);
                RHO[a] += m[b]*wab;
            }
        }
    }
    return RHO;
}



void simulation() {

    vvd R(N, vd(2));
    vvd V(N, vd(2));
    vvd VS(N, vd(2));
    vd RHO(N);
    vd S(N);
    vd m(N);
    vd C(N);


    //Initial conditions
    for (int a = 0; a < N; ++a) {
        R[a][0] = db(a)/N;
        RHO[a] = rho0;
        S[a] = s0;
        m[a] = 0.001;
    }

    // Calculate C
    for (int a = 0;  a < N; ++a) {
        C[a] = RHO[a];

        for (int b = 0; b < N; ++b) {

            if (b != a) {
                vd rab = vectordif(R[a], R[b]);
                db x = norm(rab);
                db wab = w(x);
                
                C[a] += m[b]*wab;
            }
        }
    }

    vvd VDOT(N, vd(2));
    vvd VSDOT(N, vd(2)); 
    vd SDOT(N);

    VDOT = vdot(m, RHO, R, S, V, VS);
    VSDOT = vsdot(m, RHO, R, S, V, VS);

    ofstream outFile("data.csv");
    if (!outFile) {
        cerr << "ERROR: Could not open the file. \n";
    }
    outFile << "iteration,particle,Rx,Ry,Vx,Vy,VSx,VSy,RHO\n"; // Header

    for (size_t a = 0; a < N; ++a) {
        outFile << -1 << "," << a << "," << R[a][0] << "," 
            << R[a][1] << "," << V[a][0] << "," << V[a][1] 
            << "," << VS[a][0] << "," << VS[a][1] << "," 
            << RHO[a] << "\n";
    }

    for (int i = 0; i < n; ++i) {
        if (i%5 == 0) cout << "Completed iteration " << i  << endl;

        V = matrixsum(V, matrixscale(VDOT, dt/2));
        VS = matrixsum(VS, matrixscale(VSDOT, dt/2));

        R = matrixsum(R, matrixscale(V, dt/2));

        RHO = rho(m, R, C);

        SDOT = sdot(m, RHO, R, S, V, VS);

        S = vectorsum(S, vectorscale(SDOT, dt/2));

        R = matrixsum(R, matrixscale(V, dt/2));

        RHO = rho(m, R, C);

//      We add this, we believe the paper is incomplete
        SDOT = sdot(m, RHO, R, S, V, VS);

        S = vectorsum(S, vectorscale(SDOT, dt/2));

        VDOT = vdot(m, RHO, R, S, V, VS);
        VSDOT = vsdot(m, RHO, R, S, V, VS);

        V = matrixsum(V, matrixscale(VDOT, dt/2));
        VS = matrixsum(VS, matrixscale(VSDOT, dt/2));

        // Save data

        for (size_t a = 0; a < N; ++a) {
            outFile << i << "," << a << "," << R[a][0] << "," 
                << R[a][1] << "," << V[a][0] << "," << V[a][1] 
                << "," << VS[a][0] << "," << VS[a][1] << "," 
                << RHO[a] << "\n";
        }
    }

    outFile.close();
    cout << "Data saved in 'data.csv'.\n";


}

signed main(){
    simulation();
}
