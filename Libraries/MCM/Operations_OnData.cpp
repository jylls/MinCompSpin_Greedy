#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <map>
//#include <cstring>

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
// number of binary (spin) variables:
// const unsigned int n = 9;  

// INPUT DATA FILES (optional):  
// the input datafile can also be specified directly in the main() function, as an argument of the function "read_datafile()":
// const string datafilename = "INPUT/SCOTUS_n9_N895_Data.dat"; //"INPUT/sampled.dat";

/******************************************************************************/
/***************************   Constant variables   ***************************/
/******************************************************************************/
const __int128_t un = 1;

/******************************************************************************/
/******************   TOOL Functions from "tools.cpp"   ***********************/
/******************************************************************************/
string int_to_bstring(__int128_t bool_nb, unsigned int n);
unsigned int Bitset_count(__int128_t bool_nb);

/******************************************************************************/
/***********************     READ DATA FILE    ********************************/
/******************************************************************************/
/**************    READ DATA and STORE them in Nset    ************************/
//map<__int128_t, unsigned int> read_datafile(unsigned int *N, string file = datafilename, unsigned int r=n)
map<__int128_t, unsigned int> read_datafile(unsigned int *N, string file, unsigned int r)    // O(N)  where N = data set size
{
    string line, line2;     char c = '1';
    __int128_t nb = 0, Op;
    (*N) = 0;            // N = dataset sizes
    //cout << endl << "--->> Read \"" << datafilename << "\",\t Build Nset...";

// ***** data are store in Nset:  ********************************
    map<__int128_t, unsigned int> Nset; // Nset[mu] = #of time state mu appears in the data set

    ifstream myfile (file.c_str());
    if (myfile.is_open())
    {
        while ( getline (myfile,line))
        {
            line2 = line.substr (0,r);          //take the r first characters of line
            Op = un << (r - 1);
            nb = 0;
            for (auto &elem: line2)     //convert string line2 into a binary integer
            {
                if (elem == c) { nb += Op; }
                Op = Op >> 1;
            }
            Nset[nb] += 1;
            //cout << line << endl;   //cout << nb << " :  " << int_to_bstring(nb, r) << endl;
            (*N)++;
        }
        myfile.close();
    }
    else
    {
        cout << endl << "                     ########## Unable to open file ##########" << endl << endl;
    }
    //cout << "\t\t data size N = " << (*N) << endl;
    return Nset;
}

/******************************************************************************/
/**************************     PRINT Nset   **********************************/
/******************************************************************************/
void Print_File_Nset(map<__int128_t, unsigned int> Nset, unsigned int N, unsigned int r, string OUTPUTfilename)
// map.second = nb of time that the state map.first appears in the data set
{
  map<__int128_t, unsigned int>::iterator it;
  int Ncontrol = 0;
  __int128_t un = 1;

  fstream file(OUTPUTfilename.c_str(), ios::out);
  file << "#N = " << N << endl;
  file << "#Total number of accessible states = 2^r-1 = 2^(" << r << ") - 1" << endl;
  file << "#Number of visited states, Nset.size() = " << Nset.size() << endl;
  file << "#" << endl;
  file << "#1: state \t #2: nb of pts in state \t #3: Pba state" << endl;

  for (it = Nset.begin(); it!=Nset.end(); ++it)
  {
    file << int_to_bstring((*it).first, r) << " => " << (*it).second; // << endl;
    file << "  \t  P = " << ((*it).second) / (float) N << endl;
    Ncontrol += (*it).second;
  }

  if (Ncontrol != N) { cout << "Error function \'read_Nset\': Ncontrol != N" << endl;  }

  file.close();
}

/******************************************************************************/
/*********************     CHANGE of BASIS: one datapoint  ********************/
/******************************************************************************/
// Given a choice of a model (defined by the m basis vector) --> return the new m-state (state in the new m-basis)
// Rem: must have m <= n 
__int128_t transform_mu_basis(__int128_t mu, list<__int128_t> basis)
{
  __int128_t un_i = 1, proj;
  __int128_t final_mu = 0;

  list<__int128_t>::iterator phi_i;

  for(phi_i = basis.begin(); phi_i != basis.end(); ++phi_i)
  {
    proj = (*phi_i) & mu;
    /*
    bitset<n> hi{ static_cast<unsigned long long>(proj >> 64) },
            lo{ static_cast<unsigned long long>(proj) },
            bits{ (hi << 64) | lo };

    if ( (bits.count() % 2) == 1)
    */
    if ( (Bitset_count(proj) % 2) == 1) // odd number of 1, i.e. sig_i = 1
    {
      final_mu += un_i;
    }
    un_i = (un_i << 1);
  }

  return final_mu;
}

/******************************************************************************/
/******************************   K_SET   *************************************/
/******************************************************************************/
// Build Kset for the states written in the basis of the m-chosen independent 
// operator on which the SC model is based:

map<__int128_t, unsigned int> build_Kset(map<__int128_t, unsigned int> Nset, list<__int128_t> Basis) //, bool print_bool=false)
// sig_m = sig in the new basis and cut on the m first spins 
// Kset[sig_m] = #of time state mu_m appears in the data set
{
  map<__int128_t, unsigned int>::iterator it;
  map<__int128_t, unsigned int > Kset;

  __int128_t s;        // initial state
  __int128_t sig_m;    // transformed state and to the m first spins

  unsigned int ks=0; // number of time state s appear in the dataset

  cout << endl << "--->> Build Kset..." << endl;

//Build Kset:
  for (it = Nset.begin(); it!=Nset.end(); ++it)
  {
      s = it->first;       // state s
      ks = it->second;    // # of times s appears in the data set
      sig_m = transform_mu_basis(s, Basis);
      //if (print_bool)  {  cout << int_to_bstring(s, n) << " \t" << ": \t" << int_to_bstring(sig_m, n) << endl; }

      Kset[sig_m] += ks;
  }
  cout << endl;

  return Kset;
}

/******************************************************************************/
/****************************   REDUCE K_SET   ********************************/
/******************************************************************************/
// Remove all the states that occur less than a chosen number K of times

void Reduce_Kset(map<__int128_t, unsigned int> &Kset, unsigned int K, unsigned int *N)
{
    cout << endl << "States removed from Kset:" << endl;

    map<__int128_t, unsigned int>::iterator it;
    unsigned int* counter = (unsigned int*)malloc((K+1)*sizeof(unsigned int));
    for(int i=0; i<=K; i++)
        {   counter[i]=0;   }

    for (it = Kset.begin(); it!=Kset.end(); )
    {
        if ((*it).second <= K)
        {    
            counter[(*it).second]++;
            it = Kset.erase(it);
        }
        else { ++it; }
    }

    for(int i=1; i<=K; i++)
    {   
        cout << "\t -- " << counter[i] << " states appearing ks = " << i << " times;" << endl;   
        (*N) -= i*counter[i];   // reduced the total number of states by the number of states removed
    }
    cout << endl;
}

