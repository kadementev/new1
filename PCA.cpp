#include "PCA.h"
Matrix PCA::center(Matrix m){
    double average;
    Matrix new_m(m.rows(), m.columns());
    for(int j= 0; j<m.columns(); j++){
        average= 0;
        for(int i= 0; i<m.rows(); i++){
            average+= m.get(i, j);
        }
        for(int i= 0; i<m.rows(); i++){

            new_m.set(m.get(i, j)-average/(m.rows()), i, j);
        }
    }
    return new_m;
}

Matrix PCA::autoscale(Matrix m){
    double sum;
    Matrix new_m= center(m);
    double rows= new_m.rows();
    for(int j= 0; j<new_m.columns(); j++){
        sum= 0;
        for(int i= 0; i<new_m.rows(); i++){
                sum+= new_m.get(i, j)*new_m.get(i, j);
        }
        for(int i= 0; i<m.rows(); i++){
                new_m.set(new_m.get(i, j)/sqrt((sum/(rows-1))), i, j);
        }
    }
    return new_m;
}
Matrix* PCA::NIPALS(Matrix m, int PC){
    double eps= pow(10, -8);
    Matrix D= m;
    Matrix E;
    Matrix P;
    Matrix T;
    Matrix t;
    Matrix t_old;
    Matrix p;
    Matrix d;
    Matrix* res= new Matrix[3];
    D= center(D);
    D= autoscale(D);
    E= D;
    for(int h= 0; h<PC; h++){
        t= E.slice(0,h,E.rows()-1, h);
        do{
            p= ((t.T()*E)/(t.T()*t)).T();
            p= p/p.norm();
            t_old= t;
            t= (E*p)/(p.T()*p);
            d= t_old- t;
        }while(d.norm()>eps);

        E= E- t*p.T();
        P.add_matrix(p);
        T.add_matrix(t);
    }
    res[0]=T;
    res[1]=P;
    res[2]=E;
    return res;
}
Matrix PCA::leverage(Matrix T){
    Matrix res(1, T.columns());
    Matrix temp;
    temp= T*(T.T()*T).inv()*T.T();
    for(int i=0; i<T.columns(); i++){
        res.set(temp.get(i, i),0,i);
    }
    return res;
}
Matrix PCA::deviation(Matrix E){
    Matrix res(1, E.rows());
    Matrix temp;
    for(int i= 0; i<E.rows();i++){
        temp= E.slice(i,0,i,E.columns()-1);
        temp= temp.admult(temp);
        res.set(temp.sum(),0, i);
    }
    return res;
}
double PCA::TRV(Matrix E){
    double res=0;
    Matrix temp= deviation(E);
    for(int i=0; i<E.rows(); i++){
        res+=temp.get(0, i);
    }
    return res/(E.rows()*E.columns());
}
double PCA::ERV(Matrix m, Matrix E){
    double res;
    Matrix D= center(m);
    D= autoscale(D);
    res= 1-TRV(E)/TRV(D);
    return res;
}
