#include <bitset>
#include <map>

/******************************************************************************/
/***************************   Constant variables   ***************************/
/******************************************************************************/
const unsigned int n_max = 128;  // for bitset
const unsigned int un = 1;

/******************************************************************************/
/*******************   Convert Integer to Binary string   *********************/
/******************************************************************************/
std::string int_to_bstring(__int128_t bool_nb, unsigned int n)
{
    std::string s;
    do
    {
        s.push_back( ((bool_nb & un)?'1':'0') );
    } while(bool_nb >>= 1);

    reverse(s.begin(), s.end());
    s = (std::string(n - s.length(), '0')).append(s);

    return s;
}

// Using bitset: Very bad performance (~ double the time)
/*std::string int_to_bstring_Bitset(__int128_t bool_nb)     // Very bad performance (~ double the time)
{
    std::bitset<n_max> hi{ static_cast<unsigned long long>(bool_nb >> 64) },
            lo{ static_cast<unsigned long long>(bool_nb) },
            bits{ (hi << 64) | lo };
    return bits.to_string();
}*/

/******************************************************************************/
/****************   Count number of set bits of an integer  *******************/
/******************************************************************************/
unsigned int Bitset_count(__int128_t bool_nb)     
// with Bitset: fixed time --> faster then "count_SetBits" for large values of n; for small values of n it is also faster than "count_SetBits", but not if we call systematically bitset<128>() -- as this could be a way to avoid defining "n" as a constant
// As we are more interested in optimizing for large number of variables, we will choose the option of using "bitset<128>()", instead of a manual count of the bit.
{
  std::bitset<n_max> hi{ static_cast<unsigned long long>(bool_nb >> 64) },
            lo{ static_cast<unsigned long long>(bool_nb) },
            bits{ (hi << 64) | lo };
  return bits.count();
}

unsigned int count_SetBits(__int128_t bool_nb)     // time depends on the position of the highest bit
{
    unsigned int bit_count = 0;

    while(bool_nb!=0)
    {
        if(bool_nb & un)   
            {   bit_count++;    }
        bool_nb >>= 1;
    }   
    return bit_count;
}

/******************************************************************************/
/*****************************    Check MCM   *********************************/
/******************************************************************************/
//check if *Partition* is an actual partition of the basis elements, 
// i.e., that no basis element appears in more than 1 part of the partition.
// i.e., that each basis element only appears in a single part of the partition.

std::pair<bool, unsigned int> check_partition(std::map<unsigned int, __int128_t> Partition)
{
  std::map<unsigned int, __int128_t>::iterator Part;
  __int128_t sum = 0;
  unsigned int rank = 0; 

  for (Part = Partition.begin(); Part != Partition.end(); Part++)
  {
    sum |= (*Part).second;
    rank += Bitset_count((*Part).second); 
  }

  return std::make_pair((Bitset_count(sum) == rank), rank); 
}


