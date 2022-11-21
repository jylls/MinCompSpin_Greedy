#define _USE_MATH_DEFINES
#include <iostream>
#include <sstream>
#include <fstream>
#include <bitset>
#include <map>
#include <vector>

using namespace std;

/******************************************************************************/
/*******************************    CONSTANTS     *****************************/
/******************************************************************************/
#include "data.h"

/******************************************************************************/
/***********************    READ an MCM from a FILE  **************************/
/******************************************************************************/
map<unsigned int, __int128_t> read_MCM_fromfile(string Input_MCM_file = communityfile)
{
    map<unsigned int, __int128_t> Partition;

    string line, line2;
    __int128_t Op = 1;
    Op <<= n - 1;
    vector<int> comm;

    ifstream myfile(Input_MCM_file.c_str());
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            stringstream ss(line);
            while (getline(ss, line2, '\t'))
            {
                comm.push_back(stoi(line2));
            }
            Partition[comm[1]] += Op;
            Op >>= 1;

            comm.clear();
        }
        myfile.close();
    }
    return Partition;
}

/******************************************************************************/
/*****************************    PRINT MCM   *********************************/
/******************************************************************************/
void Print_MCM_Partition(map<unsigned int, __int128_t> partition)
{
    map<unsigned int, __int128_t>::iterator it;

    for (it = partition.begin(); it != partition.end(); it++)
    {
        bitset<n> hi{ static_cast<unsigned long long>((*it).second >> 64) },
            lo{ static_cast<unsigned long long>((*it).second) },
            bits{ (hi << 64) | lo };
        cout << (*it).first << "\t " << bits << endl;
    }
    cout << endl;
}

/******************************************************************************/
/*****************************    Check MCM   *********************************/
/******************************************************************************/
//check if *Partition* is an actual partition of the r elements, 
// i.e., that no basis element appears in more than 1 part of the partition.

unsigned int count_bits(__int128_t bool_nb)
{
  bitset<n> hi{ static_cast<unsigned long long>(bool_nb >> 64) },
            lo{ static_cast<unsigned long long>(bool_nb) },
            bits{ (hi << 64) | lo };
  return bits.count();
}

pair<bool, unsigned int> check_partition(map<unsigned int, __int128_t> Partition)
{
  map<unsigned int, __int128_t>::iterator Part;
  __int128_t sum = 0;
  unsigned int rank = 0; 

  for (Part = Partition.begin(); Part != Partition.end(); Part++)
  {
    sum |= (*Part).second;
    rank += count_bits((*Part).second); 
  }

  return make_pair((count_bits(sum) == rank), rank); 
}
