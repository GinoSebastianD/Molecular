#include "iostream"
#include "algorithm"
#include "vector"
#include <string>
#include <iomanip>
#include <climits>
#include <cctype>
#include <fstream>
#include <sstream>
#include "set"
#include <functional>
using namespace std;

static const int MATCH = 1;
static const int MISMATCH = -1;
static const int GAP = -2;

inline int scoreChar(char a, char b) {
    if (a == '-' || b == '-') return GAP;
    return (a == b) ? MATCH : MISMATCH;
}

struct NWResult {
    int score;
    vector<pair<string, string>> alignments;
};

NWResult needlemanWunsch(const string& s1, const string& s2, bool enumerateAll = true) {
    int n = (int)s1.size(), m = (int)s2.size();
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));

    for (int i = 0; i <= n; i++) dp[i][0] = i * GAP;
    for (int j = 0; j <= m; j++) dp[0][j] = j * GAP;

    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++) {
            int diag = dp[i - 1][j - 1] + scoreChar(s1[i - 1], s2[j - 1]);
            int up = dp[i - 1][j] + GAP;
            int left = dp[i][j - 1] + GAP;
            dp[i][j] = max({ diag, up, left });
        }

    NWResult result;
    result.score = dp[n][m];

    set<pair<string, string>> unique_alignments;

    function<void(int, int, string, string)> traceback =
        [&](int i, int j, string a1, string a2) {
        if (i == 0 && j == 0) {
            string ra1(a1.rbegin(), a1.rend());
            string ra2(a2.rbegin(), a2.rend());
            unique_alignments.insert({ ra1, ra2 });
            return;
        }
        bool branched = false;
        if (i > 0 && j > 0 && dp[i][j] == dp[i - 1][j - 1] + scoreChar(s1[i - 1], s2[j - 1])) {
            traceback(i - 1, j - 1, a1 + s1[i - 1], a2 + s2[j - 1]);
            branched = true;
            if (!enumerateAll) return;
        }
        if (i > 0 && dp[i][j] == dp[i - 1][j] + GAP) {
            traceback(i - 1, j, a1 + s1[i - 1], a2 + '-');
            branched = true;
            if (!enumerateAll) return;
        }
        if (j > 0 && dp[i][j] == dp[i][j - 1] + GAP) {
            traceback(i, j - 1, a1 + '-', a2 + s2[j - 1]);
            branched = true;
            if (!enumerateAll) return;
        }
        (void)branched;
        };

    traceback(n, m, "", "");
    result.alignments.assign(unique_alignments.begin(), unique_alignments.end());
    return result;
}

struct StarMSAResult {
    int centerIndex;
    vector<string> alignedSequences;
    int centerStarScore;
    long long spScore;
};

void mergeIntoMSA(vector<string>& msaRows, vector<int>& mergedIdx, int centerRow,
    const string& alignedCenterNew, const string& alignedOtherNew, int otherIndex) {
    const string& alignedCenterOld = msaRows[centerRow];

    vector<string> newRows(msaRows.size());
    string newOtherRow;

    size_t p1 = 0, p2 = 0;
    while (p1 < alignedCenterOld.size() || p2 < alignedCenterNew.size()) {
        char c1 = (p1 < alignedCenterOld.size()) ? alignedCenterOld[p1] : '\0';
        char c2 = (p2 < alignedCenterNew.size()) ? alignedCenterNew[p2] : '\0';

        if (p1 < alignedCenterOld.size() && c1 == '-') {
            for (size_t r = 0; r < msaRows.size(); r++) newRows[r] += msaRows[r][p1];
            newOtherRow += '-';
            p1++;
        }
        else if (p2 < alignedCenterNew.size() && c2 == '-') {
            for (size_t r = 0; r < msaRows.size(); r++) newRows[r] += '-';
            newOtherRow += alignedOtherNew[p2];
            p2++;
        }
        else {
            for (size_t r = 0; r < msaRows.size(); r++) newRows[r] += msaRows[r][p1];
            newOtherRow += alignedOtherNew[p2];
            p1++; p2++;
        }
    }

    msaRows = newRows;
    msaRows.push_back(newOtherRow);
    mergedIdx.push_back(otherIndex);
}

StarMSAResult starMSA(const vector<string>& seqs, bool enumerateAllPairwise = false) {
    int k = (int)seqs.size();
    vector<vector<int>> pairScore(k, vector<int>(k, 0));
    vector<vector<pair<string, string>>> pairAlign(k, vector<pair<string, string>>(k));
    vector<vector<int>> optimalCount(k, vector<int>(k, 0));

    for (int i = 0; i < k; i++)
        for (int j = 0; j < k; j++) {
            if (i == j) continue;
            NWResult r = needlemanWunsch(seqs[i], seqs[j], enumerateAllPairwise);
            pairScore[i][j] = r.score;
            optimalCount[i][j] = (int)r.alignments.size();
            pairAlign[i][j] = r.alignments[0];
        }

    if (enumerateAllPairwise) {
        cout << "  Numero de alineamientos optimos por par (empates en NW):\n";
        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++)
                cout << "    S" << (i + 1) << " - S" << (j + 1)
                << " : score=" << pairScore[i][j]
                << "  (" << optimalCount[i][j] << " alineamiento(s) optimo(s))\n";
    }

    int center = 0, bestSum = INT_MIN;
    for (int i = 0; i < k; i++) {
        int sum = 0;
        for (int j = 0; j < k; j++) if (i != j) sum += pairScore[i][j];
        if (sum > bestSum) { bestSum = sum; center = i; }
    }

    vector<string> msaRows = { seqs[center] };
    vector<int> mergedIdx = { center };

    for (int j = 0; j < k; j++) {
        if (j == center) continue;
        mergeIntoMSA(msaRows, mergedIdx, 0, pairAlign[center][j].first, pairAlign[center][j].second, j);
    }

    vector<string> aligned(k);
    for (size_t r = 0; r < msaRows.size(); r++) aligned[mergedIdx[r]] = msaRows[r];

    long long sp = 0;
    int L = (int)aligned[0].size();
    for (int c = 0; c < L; c++)
        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++)
                sp += scoreChar(aligned[i][c], aligned[j][c]);

    StarMSAResult res;
    res.centerIndex = center;
    res.alignedSequences = aligned;
    res.centerStarScore = bestSum;
    res.spScore = sp;
    return res;
}

void printMSA(const vector<string>& names, const vector<string>& aligned) {
    size_t maxName = 0;
    for (auto& n : names) maxName = max(maxName, n.size());
    for (size_t i = 0; i < aligned.size(); i++) {
        cout << "  " << left << setw((int)maxName + 2) << names[i] << aligned[i] << "\n";
    }
}

void ejemploClase() {
    cout << "\n";
    cout << " EJEMPLO 1: MSA Estrella sobre las 5 cadenas vistas en clase\n";
    cout << "\n";
    vector<string> names = { "S1", "S2", "S3", "S4", "S5" };
    vector<string> seqs = {
        "ATTGCCATT",
        "ATGGCCATT",
        "ATCCAATTTT",
        "ATCTTCTT",
        "ACTGACC"
    };

    cout << "Secuencias de entrada:\n";
    printMSA(names, seqs);
    cout << "\n";

    StarMSAResult res = starMSA(seqs, true);

    cout << "\nSecuencia centro elegida: S" << (res.centerIndex + 1)
        << " (mayor suma de similitud con las demas)\n\n";

    cout << "Alineamiento multiple resultante (metodo Estrella):\n";
    printMSA(names, res.alignedSequences);

    cout << "\nScore 'centro-estrella' (suma NW(centro, otros)): " << res.centerStarScore << "\n";
    cout << "Score suma-de-pares (SP) del MSA final           : " << res.spScore << "\n";
    cout << "(match=" << MATCH << ", mismatch=" << MISMATCH << ", gap=" << GAP << ")\n\n";
}

bool leerSecuencias(const string& ruta, vector<string>& nombres, vector<string>& secuencias) {
    ifstream archivo(ruta);
    if (!archivo) return false;

    string linea;
    string muestraActual;

    auto extraerBases = [](const string& texto) {
        string bases;
        for (char caracter : texto) {
            char base = static_cast<char>(toupper(static_cast<unsigned char>(caracter)));
            if (base == 'A' || base == 'C' || base == 'G' || base == 'T' || base == 'N') {
                bases += base;
            }
        }
        return bases;
    };

    while (getline(archivo, linea)) {
        if (!linea.empty() && linea.back() == '\r') linea.pop_back();
        if (linea.empty()) continue;

        size_t marcadorForward = linea.find("F:");
        size_t marcadorReverse = linea.find("R:");
        size_t marcador = marcadorForward != string::npos ? marcadorForward : marcadorReverse;

        if (marcador != string::npos) {
            if (marcadorForward != string::npos) {
                stringstream entrada(linea.substr(0, marcador));
                entrada >> muestraActual;
            }

            string bases = extraerBases(linea.substr(marcador + 2));
            if (!bases.empty() && !muestraActual.empty()) {
                nombres.push_back(muestraActual + (marcadorForward != string::npos ? "_F" : "_R"));
                secuencias.push_back(bases);
            }
        }
    }

    return !secuencias.empty();
}

void ejemploBRCA1() {
    cout << "\n";
    cout << " EJEMPLO 2: MSA Estrella sobre muestras BRCA1 (BRCA1.txt)\n";
    cout << "\n";

    string ruta;
    vector<string> nombres, secuencias;
    cout << "Ruta del archivo: ";
    cin >> ruta;

    if (!leerSecuencias(ruta, nombres, secuencias)) {
        cout << "No se pudo leer el archivo o no contiene secuencias.\n";
        return;
    }

    vector<string> namesF, seqsF, namesR, seqsR;
    for (size_t i = 0; i < secuencias.size(); i++) {
        if (nombres[i].find("_R") != string::npos) {
            namesR.push_back(nombres[i]);
            seqsR.push_back(secuencias[i]);
        }
        else {
            namesF.push_back(nombres[i]);
            seqsF.push_back(secuencias[i]);
        }
    }

    if (namesR.empty()) {
        namesR = namesF;
        seqsR = seqsF;
    }

    cout << "--- Cadenas FORWARD ---\n";
    cout << "Secuencias de entrada:\n";
    printMSA(namesF, seqsF);
    StarMSAResult resF = starMSA(seqsF, false);
    cout << "\nSecuencia centro elegida: " << namesF[resF.centerIndex] << "\n\n";
    cout << "Alineamiento multiple resultante:\n";
    printMSA(namesF, resF.alignedSequences);
    cout << "\nScore centro-estrella: " << resF.centerStarScore
        << "   |   Score SP: " << resF.spScore << "\n\n";

    cout << "--- Cadenas REVERSE ---\n";
    cout << "Secuencias de entrada:\n";
    printMSA(namesR, seqsR);
    StarMSAResult resR = starMSA(seqsR, false);
    cout << "\nSecuencia centro elegida: " << namesR[resR.centerIndex] << "\n\n";
    cout << "Alineamiento multiple resultante:\n";
    printMSA(namesR, resR.alignedSequences);
    cout << "\nScore centro-estrella: " << resR.centerStarScore
        << "   |   Score SP: " << resR.spScore << "\n\n";

    cout << "Interpretacion: al alinear solo forward entre si (y solo reverse entre si),\n"
        << "se comparan primers que amplifican DISTINTOS exones/regiones del gen BRCA1;\n"
        << "por eso los scores tienden a ser bajos y el alineamiento muestra pocos matches\n"
        << "reales y muchos gaps/mismatches: no comparten homologia de secuencia significativa,\n"
        << "solo el hecho de pertenecer al mismo gen y protocolo de amplificacion LR-PCR.\n";
}

int main() {
    int opcion;

    do {
        cout << "\n MENU \n";
        cout << "1. Ejecutar ejemplo 1\n";
        cout << "2. Ejecutar ejemplo 2\n";
        cout << "0. Salir\n";
        cout << "Seleccione una opcion: ";
        cin >> opcion;

        switch (opcion) {
        case 1:
            ejemploClase();
            break;
        case 2:
            ejemploBRCA1();
            break;
        case 0:
            cout << "Programa finalizado.\n";
            break;
        default:
            cout << "Opcion no valida.\n";
        }
    } while (opcion != 0);

    return 0;
}