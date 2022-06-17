#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include "Matrix.h"
#include "PCA.h"
using namespace std;
int main() {
    Matrix* res;
    Matrix data;
    PCA pca;
    int PC= 4; //1<PC<13
    ifstream file("C:\\data.txt");
    data<<file;
    res= pca.NIPALS(data, PC);
    cout<<"PCA(PC= "<<PC<<"):"<<endl;
    cout<<"Scores:"<<endl<<res[0]<<"Lodings:"<<endl<<res[1];
    cout<<"Leverage:"<<endl<<pca.leverage(res[0])<<"Deviation"<<endl<<pca.deviation(res[2]);
    cout<<"TRV: "<<pca.TRV(res[2])<<" ERV: "<<pca.ERV(data,res[2])<<endl;
}
