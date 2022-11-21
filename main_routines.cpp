#include <map>
#include <list>
#include <fstream>
#include <sstream>
#include <ctime> // for chrono
#include <ratio> // for chrono
#include <chrono> // for chrono

using namespace std;
using namespace std::chrono;

/******************************************************************************/
/**********************    CONSTANTS AND FUNCTIONS    *************************/
/******************************************************************************/
#include "data.h"
#include "library.h"

/******************************************************************************/
/**********************   ROUTINES for GREEDY SEARCH   ************************/
/******************************************************************************/
map<unsigned int, __int128_t> MCM_GreedySearch(map<__int128_t, unsigned int> Kset, unsigned int N, unsigned int r = n);

map<unsigned int, __int128_t> MCM_GreedySearch_AND_printInfo(map<__int128_t, unsigned int> Kset, unsigned int N, unsigned int r = n)
{
    cout << "######### START GREEDY SEARCH #########" << endl;
    // *** Calculate the optimal partition
    auto start = chrono::system_clock::now();
    map<unsigned int, __int128_t> fp1 = MCM_GreedySearch(Kset, N, r);
    auto end = chrono::system_clock::now();

    // *** Time it takes to find partition
    chrono::duration<double> elapsed = end - start;

    cout << "######### EMPERICAL #########" << endl;
    // Entropy of dataset
    double H = Entropy(Kset, N);
    cout << "H : " << H << ". Range: [0, " << n << "]" << endl << endl;

    cout << "#########  GREEDY   #########" << endl;
    // Log evidence of MCM
    double LE_g = LogE_MCM(Kset, fp1, N);
    Print_MCM_Partition(fp1);

    cout << "Elapsed time      : " << elapsed.count() << "s" << endl;
    cout << "Log-evidence      : " << LE_g << endl;
    cout << "Average comm size : " << (double)n / (double)fp1.size() << endl << endl;

    return fp1;
}

/******************************************************************************/
/********************   GREEDY SEARCH in ORIGINAL BASIS   *********************/
/******************************************************************************/
map<unsigned int, __int128_t> MCM_GreedySearch_OriginalBasis(map<__int128_t, unsigned int> Nset, unsigned int N)
{
    return MCM_GreedySearch(Nset, N, n);
}

map<unsigned int, __int128_t> MCM_GreedySearch_ChosenBasis(map<__int128_t, unsigned int> Nset, unsigned int N, list<__int128_t> Basis_li)
{
    map<__int128_t, unsigned int> Kset = build_Kset(Nset, Basis_li, false);
    return MCM_GreedySearch(Kset, N, Basis_li.size());
}

/******************************************************************************/
/********************   GREEDY SEARCH for REDUCED DATASET  ********************/
/******************************************************************************/
map<unsigned int, __int128_t> MCM_ReducedGreedySearch_AND_PrintInfo(map<__int128_t, unsigned int> Kset, unsigned int K, unsigned int N, unsigned int r = n)
{
    unsigned int N_reduced = N;
    map<__int128_t, unsigned int> Kset_reduced = Kset;

    cout << "######### REDUCE Kset TO STATES OCCURING AT LEAST K TIMES: K = " << K << "  #########" << endl;
    Reduce_Kset(Kset_reduced, K, &N_reduced);

    cout << "After reduction: " << endl;
    cout << "\t Number of datapoints = " << N_reduced << endl;
    cout << "\t Number of different observed states = " << Kset_reduced.size() << endl;

    cout << "######### START GREEDY SEARCH #########" << endl;
    // *** Calculate the optimal partition
    auto start = chrono::system_clock::now();
    map<unsigned int, __int128_t> fp_reduced = MCM_GreedySearch(Kset_reduced, N_reduced, r);
    auto end = chrono::system_clock::now();

    // *** Time it takes to find partition
    chrono::duration<double> elapsed = end - start;

    cout << "######### EMPERICAL #########" << endl;
    // Entropy of dataset
    double H = Entropy(Kset_reduced, N_reduced);
    cout << "Entropy, H : " << H << ". Range: [0, " << n << "]" << endl << endl;

    cout << "#########  MCM GREEDY INFO  #########" << endl;
    // Log evidence of MCM
    double LE_reduced = LogE_MCM(Kset_reduced, fp_reduced, N_reduced);
    double LE = LogE_MCM(Kset, fp_reduced, N);
    Print_MCM_Partition(fp_reduced);

    cout << "Elapsed time      : " << elapsed.count() << "s" << endl;
    cout << "Log-evidence reduced     : " << LE_reduced << endl;
    cout << "Log-evidence original data     : " << LE << endl;
    cout << "Average comm size : " << (double)n / (double)fp_reduced.size() << endl << endl;

    Kset_reduced.clear();

    return fp_reduced;
}

/******************************************************************************/
/*************************   READING MCM from a FILE   ************************/
/******************************************************************************/
map<unsigned int, __int128_t> read_MCM_fromfile(string Input_MCM_file = communityfile);

map<unsigned int, __int128_t> read_MCM_fromfile_AND_printInfo(map<__int128_t, unsigned int> Kset, unsigned int N, string Input_MCM_file = communityfile)
{
    cout << "#########  THEORETICAL   #########" << endl;
    map<unsigned int, __int128_t> fp2 = read_MCM_fromfile(Input_MCM_file);

    double LE_t = LogE_MCM(Kset, fp2, N);
    Print_MCM_Partition(fp2);

    cout << "Log-evidence      : " << LE_t << endl;
    cout << "Average comm size : " << (double)n / (double)fp2.size() << endl << endl;

    return fp2;
}

/******************************************************************************/
/****************************   COMPARING TWO MCMs   **************************/
/******************************************************************************/
void compare_two_MCMs_AND_printInfo(map<__int128_t, unsigned int> Kset, unsigned int N, map<unsigned int, __int128_t> fp1, map<unsigned int, __int128_t> fp2)
{
    cout << "#########  COMPARATIVE MEASURES   #########" << endl;
    double VOI = Var_of_Inf(fp1, fp2);
    double NMI = Norm_Mut_info(fp1, fp2);
    string istrue = is_subset(fp1, fp2) ? "Yes" : "No";

    double LE_1 = LogE_MCM(Kset, fp1, N);
    double LE_2 = LogE_MCM(Kset, fp2, N);

    cout << "Is MCM_1 \'subset\' of MCM_2    : " << istrue << endl;
    cout << "Variation of Information      : " << VOI << endl;
    cout << "Normalized Mutual Information : " << NMI << endl;
    cout << "Difference in Log-Evidence    : LogE1-LogE2 = " << LE_1 - LE_2 << endl << endl;
}



