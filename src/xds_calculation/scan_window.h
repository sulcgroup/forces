#include <iostream>
#include <cstdio>
#include <string>
#include <regex>
#include <string>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <cstdlib>
#include <exception>
#include <utility>      
#include <stdexcept>
#include <algorithm>
#include <locale>
#include <cmath>

using namespace std;


typedef pair<string,string> name_seq;
//--------------------------------------------------------------------------------------------------------------------------------------
list<name_seq>  read_fasta( const string& fname );


//--------------------------------------------------------------------------------------------------------------------------------------
class AssignSequences {
    string _seq;
    int _window_size;
    int _shift_size;
    int *_matrix;
    bool _wobble; 
    int _start; 
    int _stop;
    int _N; //matrix size
   // list<double> xds;

  public:
   AssignSequences(const string& seq, int window_size, int shift_size,bool wobble=true) : _seq(seq), _window_size(window_size), _shift_size(shift_size), _wobble(wobble) 
   {
       _N = _window_size;
       this->_matrix = nullptr; // new int [N*N] ();
   }
   
   ~AssignSequences() {delete [] this->_matrix;}

  
    //helper functions
   bool are_compatible(char a, char b);
   
   int find_longest_stretch(const string& seq, int &sA_start, int& sA_end, int& sB_start, int& sB_end, bool include_wobble=true);
   static double get_xds(const string& myseq, int start, int end, int maxlen,double c = -2.2);
   //---------------------------------------------------------------------------------------------------------------------------------------  
   bool check_segments(string &A, string &B);


   void fill_matrix(int seq_start, int seq_end=-1); //seq_end = seq_start + window_size
   //int extract_ds_from_matrix(int &sA_start, int& sA_end, int& sB_start, int& sB_end);
   double extract_xds_from_matrix(int &maxlen,int& SaS, int& SaE, int& SbS, int& SbE);
   int shift_and_fill_matrix(int seq_start,int shift); //seq_start is the position where the new window starts, and it is shiftd by shift from old seq start 
   int get_maxlen_from_matrix(int &sA_start, int& sA_end, int& sB_start, int& sB_end);
   
   void scan_and_do_xds(ostream &out);
   void scan_and_do_maxlen(ostream &out);
   
   

} ;

void process_sequence_file(const char *fastafile,int window_size,int step_size,string outputfilename);

//------------------------------------------------------------------------------------------------------------------------------------
/*
int assign_xds_scanning_window(const string& seq, int &sA_start, int& sA_end, int& sB_start, int& sB_end, bool include_wobble=true)
{
 int N = seq.length();
 int *matrix = new  int [N*N] ();
 int maxlen = 0;
 int res_i = 0;
 int res_j = 0;
 for(int i = 0; i < N; i++)
  {
      for(int j = i; j < N; j++)
      {
          if (are_compatible(seq[i],seq[j],include_wobble))
          {
              matrix[i + j*N] = 1;
          }
      }
  }

  for(int i = 0; i < N; i++)
  {
      for(int j = i; j < N; j++)
      {
          if (matrix[i + j*N])
          {
              int sublen = 1;
              for(int k = 1; i-k >= 0 && k+j < N; k++)
              {
                  if( matrix[ (i-k) + (j+k)*N])
                  {
                      sublen += 1;
                  }
                  else
                  {
                      break;
                  }
                  
              }
              if(sublen > maxlen)
              {
                  maxlen = sublen;
                  res_i = i;
                  res_j = j;
              } 
          }
      }
  }
  
 sA_start  = res_i-maxlen+1;
 sA_end    = res_i;
 sB_start  = res_j;
 sB_end    = res_j + maxlen-1;

 delete [] matrix;
 return maxlen;
}
*/