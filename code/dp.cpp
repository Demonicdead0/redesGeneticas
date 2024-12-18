#include<bits/stdc++.h>
#include<fstream>

using namespace std;

#define db(v) for (auto x: v) cout<<x<<" "; cout<<endl;


/*
===============================
ENVIAR DATOS
==============================
*/

void guardar(vector<vector<int>>& datos, const string & nombre){
    ofstream archivo(nombre);
    if (!archivo){
        cerr << "error"<<endl;
        return;
    }

    for (const auto&fila: datos){
        for (int i = 0; i < fila.size(); i++){
            archivo << fila[i];
            if (i < fila.size() -1){
                archivo <<" ";
            }
        }
        archivo << "\n"<<endl;
    }

    archivo.close();

}
/*
===============================
INFORMACION MUTUA
==============================
 */

// Función para calcular la frecuencia de cada valor en un vector
std::map<int, double> calcularFrecuencias(const std::vector<int>& datos) {
    std::map<int, double> frecuencias;
    for (int valor : datos) {
        frecuencias[valor]++;
    }
    for (auto& par : frecuencias) {
        par.second /= datos.size();
    }
    return frecuencias;
}

// Función para calcular la entropía de un conjunto de datos
double calcularEntropia(const std::vector<int>& datos) {
    auto frecuencias = calcularFrecuencias(datos);
    double entropia = 0.0;
    for (const auto& par : frecuencias) {
        double p = par.second;
        if (p > 0) {
            entropia -= p * std::log2(p);
        }
    }
    return entropia;
}

// Función para calcular la entropía conjunta de dos variables
double calcularEntropiaConjunta(const std::vector<int>& X, const std::vector<int>& Y) {
    std::map<std::pair<int, int>, double> frecuenciasConjuntas;
    for (size_t i = 0; i < X.size(); ++i) {
        frecuenciasConjuntas[{X[i], Y[i]}]++;
    }
    for (auto& par : frecuenciasConjuntas) {
        par.second /= X.size();
    }
    double entropiaConjunta = 0.0;
    for (const auto& par : frecuenciasConjuntas) {
        double p = par.second;
        if (p > 0) {
            entropiaConjunta -= p * std::log2(p);
        }
    }
    return entropiaConjunta;
}

// Función para calcular la Información Mutua
double calcularInformacionMutua(const std::vector<int>& X, const std::vector<int>& Y) {
    double H_X = calcularEntropia(X);
    double H_Y = calcularEntropia(Y);
    double H_XY = calcularEntropiaConjunta(X, Y);
    return H_X + H_Y - H_XY;
}

// Función para calcular la Información Mutua Penalizada
double calcularInformacionMutuaPenalizada(const std::vector<int>& X, const std::vector<int>& Y) {
    double MI = calcularInformacionMutua(X, Y);
    std::set<int> valoresX(X.begin(), X.end());
    std::set<int> valoresY(Y.begin(), Y.end());
    double penalizacion = 0.5 * std::log2(valoresX.size() * valoresY.size());
    return MI - penalizacion;
}

double calcularMI(vector<vector<int>> & datos, int gen, vector<int>& genes_influentes){
    int nro_genes = datos[0].size();
    int tiempo = datos.size();
    map<vector<int>, map<int,int>> tabla;
   // db(genes_influentes);
    for (int ti = 1; ti < tiempo; ti++){
        vector<int> caso;
        for (auto gen_influente: genes_influentes){
            caso.push_back(datos[ti-1][gen_influente]);
        }
        //db(caso);
        tabla[caso][datos[ti][gen]]++;
    }
    vector<int> x,y;
    for (auto xx: tabla){
        x.push_back(xx.second[0]);
        y.push_back(xx.second[1]);
    }
    //cout<<"----"<<endl;
    //db(x);
    ///db(y);
    return calcularInformacionMutua(x,y);
}

/*
===============================
FUNCIONES BOOLENAS
==============================
 */

int predecir_gen(vector<vector<int>>& datos, int gen, vector<int> & genes_influentes, vector<int> & estado_actual){
    int nro_genes = datos[0].size();
    int tiempo = datos.size();
    map<vector<int>, map<int,int>> tabla;
    //db(genes_influentes);
    for (int ti = 1; ti < tiempo; ti++){
        vector<int> caso;
        for (auto gen_influente: genes_influentes){
            caso.push_back(datos[ti-1][gen_influente]);
        }
        //db(caso);
        tabla[caso][datos[ti][gen]]++;
    }
    vector<int> estado;
    for (auto x: genes_influentes) estado.push_back(estado_actual[x]);
    int x = tabla[estado][0];
    int y = tabla[estado][1];
    if (x == y) return rand()%2;
    if (x > y){
        return 0;
    }else{
        return 1;
    }
}

/*
===============================
BRUTE FORCE
==============================
 */

vector<vector<int>> dp(vector<vector<int>> & datos){
    cerr<<"DP"<<endl;
    vector<vector<int>> genes_caracteristicos;
    int nro_genes = datos[0].size();
    for (int i = 0; i < nro_genes; i++){
        int gen = i;
        int nro_genes = datos[0].size();
        vector<vector<double>> dp(nro_genes+1, vector<double>(nro_genes+1,1e9));
        vector<vector<vector<int>>> opt(nro_genes+1, vector<vector<int>> (nro_genes +1, vector<int>(0)));
        dp[0][0] = calcularMI(datos, gen, opt[0][0]);
        for (int ii = 0; ii < nro_genes; ii++){
            for (int jj = 0; jj < nro_genes; jj++) if (ii >= jj){
                if (dp[ii+1][jj] >= dp[ii][jj]){
                    dp[ii+1][jj] = dp[ii][jj];
                    opt[ii+1][jj] = opt[ii][jj];
                }
                vector<int> part = opt[ii][jj];
                part.push_back(ii);
                int subans = calcularMI(datos, gen, part);
                if (dp[ii+1][jj+1] >= subans){
                    dp[ii+1][jj+1] = subans;
                    opt[ii+1][jj] = part;
                }
            }
        }
        pair<double, vector<int>> mejor;
        mejor.first = 1e9;
        mejor.second = vector<int>(0);
        for (int ii = 1; ii <= nro_genes; ii++){
            if (mejor.first >= dp[nro_genes][ii]){
                mejor.first = dp[nro_genes][ii];
                mejor.second = opt[nro_genes][ii];
            }
        }
        genes_caracteristicos.push_back(mejor.second);
    }
    return genes_caracteristicos;
}




/*
===============================
MAIN
==============================
 */
int32_t main(){
    //ios::sync_with_stdio(false), cout.tie(0), cin.tie(0);
    int tiempo, nro_genes;
    double none; cin>>tiempo>>nro_genes>>none;
    cout<<tiempo<<" "<<nro_genes<<" "<<none<<endl;
    vector<vector<int>> genes_tiempo(tiempo +1, vector<int> (nro_genes));
    for (int i = 0; i <= tiempo;i++){
        for (int j = 0; j < nro_genes; j++){
            cin>>genes_tiempo[i][j];
            //cout<<i<<endl;
            //cout<<genes_tiempo[i][j]<<endl;
        }
            
    }

    
    vector<vector<int>> grafo1 = dp(genes_tiempo);
   
    
    // construir los datos 
    vector<vector<int>> nuevos_datos(tiempo+1, vector<int>(nro_genes));
    nuevos_datos[0] = genes_tiempo[0];
    for (int i = 1; i < tiempo; i++){
        for (int j = 0; j < nro_genes; j++){
            nuevos_datos[i][j] = predecir_gen(genes_tiempo, j, grafo1[j], nuevos_datos[i-1]); 
        }
    }
    for (int i = 0; i < tiempo; i++){
        for (int j = 0; j < nro_genes; j++) cout<<nuevos_datos[i][j]<<" "; 
        if (i != tiempo -1) cout<<endl;
    }
}