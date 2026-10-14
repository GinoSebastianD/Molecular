#include "iostream"
#include "vector"

#include "algorithm"
#include "string"


using namespace std;

const int igual = 1;
const int desigual = -1;
const int gap = -2;

string s, t;

vector<vector<int>> M;
vector<pair<string,string>> results;

int iguales(char a, char b){
    return (a == b) ? igual : desigual;
}


void retroceder(int m , int n, string& alS, string& alT){

    if(m == 0 && n == 0){
        return results.push_back( {string(alS.rbegin() , alS.rend()) , string(alT.rbegin(),alT.rend())} );
        return;
    }

    if( m > 0 && n > 0 && M[m][n] == M[m-1][n-1] + iguales(s[m-1],t[n-1] )){
        alS.push_back(s[m-1]);
        alT.push_back(t[n-1]);

        retroceder(m-1,n-1,alS,alT );
        alS.pop_back();
        alT.pop_back();
    }
    if( m >0 && M[m][n] == M[m-1][n] + gap){
        alS.push_back(s[m-1]);
        alT.push_back('-');
        retroceder(m-1,n,alS,alT);
        alS.pop_back();
        alT.pop_back();
    }

    if(n > 0 && M[m][n]== M[m][n-1] + gap){
        alS.push_back('-');
        alT.push_back(t[n-1]);
        retroceder(m,n-1,alS,alT);
        alS.pop_back();
        alT.pop_back();

    }



}


int main(){

    getline(cin,s);
    getline(cin,t);
    int m = s.size();
    int n = t.size();

    M.assign(m+1,vector<int>(n+1,0));

    for(int i = 1 ; i <= m ; i++){
        M[i][0] = i * gap;
    }
     for(int j = 1 ; j <= n ; j++){
        M[0][j] = j * gap;
    }

    for(int i = 1; i <= m; i++){
        for(int j = 1; j <= n; j++){
            M[i][j] = max({M[i-1][j-1] + iguales(s[i-1],t[j-1]) , M[i-1][j] + gap , M[i][j-1]+ gap}   );
        }
    }

    cout<<"\nScore = " << M[m][n] << "\n";

    string alS, alT;
    retroceder(m,n,alS,alT);
    cout<<"\n";
    cout <<results.size();




}
