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
    return calcularInformacionMutuaPenalizada(x,y);
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

vector<vector<int>> bfb(vector<vector<int>>& datos){
    cerr<<"BFB"<<endl;
    vector<vector<int>> genes_caracteristicos;
    int nro_genes = datos[0].size();
    for (int i = 0; i < nro_genes; i++){
        int gen = i;
        // start bfb
        vector<int> subopt;
        pair<double, vector<int>> mejor;
        mejor.first = 1e9;
        mejor.second = vector<int>(0);
        while(!subopt.empty()){
            int gen_seleccionado = -1;
            for (auto elemento: subopt){
                vector<int> suboptimo;
                for (int ii = 0; ii < subopt.size(); ii++) if (subopt[ii] != elemento) suboptimo.push_back(ii);
                double subans = calcularMI(datos, gen, suboptimo);
                if (mejor.first >= subans){
                    mejor.first =subans;
                    mejor.second = suboptimo;
                    gen_seleccionado = elemento;
                }
            }
            if (gen_seleccionado == -1) break;
            else{
                vector<int> suboptimo;
                for (int ii = 0; ii < subopt.size(); ii++) if (subopt[ii] != gen_seleccionado) suboptimo.push_back(ii);
                subopt = suboptimo;
            }
        }
        genes_caracteristicos.push_back(subopt);
    }
    return genes_caracteristicos;
}

vector<vector<int>> sfs(vector<vector<int>>& datos){
    cerr<<"SFS"<<endl;
    vector<vector<int>> genes_caracteristicos;
    int nro_genes = datos[0].size();
    for (int i = 0; i < nro_genes; i++){
        int gen = i;
        // start sfs
        set<int> genes_restantes;
        for (int ii = 0; ii < nro_genes; ii++) genes_restantes.insert(ii);
        vector<int> subopt;
        pair<double, vector<int>> mejor;
        mejor.first = 1e9;
        mejor.second = vector<int>(0);
        while (!genes_restantes.empty()){
            int gen_seleccionado = -1;
            for (auto elemento: genes_restantes){
                subopt.push_back(elemento);
                double subans = calcularMI(datos,gen, subopt);
                if (mejor.first >= subans){
                    mejor.first = subans;
                    mejor.second = subopt;
                    gen_seleccionado = elemento;
                }
                subopt.pop_back();
            }
            if (gen_seleccionado == -1) break;
            else{
                subopt.push_back(gen_seleccionado);
                genes_restantes.erase(genes_restantes.find(gen_seleccionado));
            }
        }

        genes_caracteristicos.push_back(subopt);
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

    
    vector<vector<int>> grafo1 = bfb(genes_tiempo);
   
    
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