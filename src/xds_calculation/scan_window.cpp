
#include "scan_window.h"
#include <algorithm>

using namespace std;

 
const std::string WHITESPACE = " \n\r\t\f\v";
 
std::string trim(const std::string s) {
    size_t pos = s.find_first_of(WHITESPACE);
    if (pos != std::string::npos) {
        return s.substr(0, pos);
    }
    return s;
}

//--------------------------------------------------------------------------------------------------------------------------------------
list<name_seq>  read_fasta( const string& fname ){
   
    std::ifstream input(fname);
    
    list<name_seq> loaded_seqs;

    if(!input.good()){
        std::cerr << "Error opening '"<< fname <<"'. Bailing out." << std::endl;
        throw runtime_error("Could not open file");
    }
 
    std::string line, name, content;
    while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Print out what we read from the last entry
                //std::cout << name << " : " << content << std::endl;
                for (auto & c: content) c = toupper(c);
                loaded_seqs.push_back( make_pair(name,content));
                name.clear();
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }
    if( !name.empty() ){ // Print out what we read from the last entry
        //std::cout << name << " : " << content << std::endl;
        //content = std::toupper(content.being(,content.end()));
        for (auto & c: content) c = toupper(c);
        loaded_seqs.push_back( make_pair(name,content));
    }

    //cout << "Loaded " << loaded_seqs.size() << " sequences " << endl;
    //loaded_seqs.push_back( make_pair(name,content));
    return loaded_seqs;
}
//---------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------
void AssignSequences::fill_matrix(int seq_start, int seq_end)
{
    if (seq_start < 0 || seq_start >= _seq.length())
    {
        throw string("Incorrect seq_start");
    }
    if(seq_end != -1)
    {
        this->_N = seq_end - seq_start + 1;
    }

    this->_start = seq_start;
    this->_stop  = seq_end;

    if(this->_matrix == nullptr)
    {
        this->_matrix = new int [_N * _N] ();
    }

    int maxlen = 0;
    int res_i = 0;
    int res_j = 0;
    for(int i = 0; i < _N; i++)
    {
        for(int j = i; j < _N; j++)
        {
            if (this->are_compatible(_seq[i+seq_start],_seq[j+seq_start]) )
            {
                _matrix[i + j*_N] = 1;
            }
            else
            {
                 _matrix[i + j*_N] = 0;
            }
            
        }
    }   
}
//-------------------------------------------------------------------------------------------------------------------------------------
int AssignSequences::shift_and_fill_matrix(int seq_start,int shift) //seq_start is the position where the new window starts, and it is shiftd by shift from old seq start 
{
 //cout << " Doing shift " <<  seq_start << " and " << shift << endl;
 //cout << "start and stop: " << this->_start << " " << this->_stop << endl;
 if (shift + this->_stop > this->_seq.length())
 {
     shift = this->_seq.length() - this->_stop - 1;
     if (shift <= 0)
     {
         return -1;
     }
 }
 
 for(int i = shift; i < _N; i++)
 {
    for(int j = i; j < _N; j++)
    {
      if(i-shift + (j-shift)*_N  >= _N * _N  || i-shift + (j-shift)*_N  < 0)
      {
           cerr << "Outside the number boundaries " << i-shift  << " " << (j-shift)  << " should be " <<  _N;
           throw string(" trying to shift outside boundaries ");

      }
      if(i+ j*_N  >= _N * _N  || i + j*_N  < 0)
      {
           cerr << "2) Outside the number boundaries " << i  << " " << j  << " should be " <<  _N;
           throw string(" trying to shift outside boundaries ");

      }
      this->_matrix[i-shift + (j-shift)*_N ] =  this->_matrix[i + j*_N]  ;
    }
 }
 
 this->_start += shift;
 if(this->_start != seq_start)
 {
     cout << "? misalignment? " << this->_start <<  " " << seq_start << endl;
 }
 this->_stop += shift;

 for(int i = 0; i < _N; i++)
 {
        for(int j = i < _N - shift ? _N - shift : i ; j < _N; j++)
        {
            //int seq_i = i;
            //int seq_j = j; //;;;;; - _N - shift - 1; 
            if(i + j*_N >= _N * _N  ||  i + j*_N  < 0)
            {
                cerr << "Outside the number boundaries " << i  << " " << j  << " should be " <<  _N;
                throw string(" trying to shift outside boundaries ");

            }
            if (this->are_compatible(_seq[i+_start],_seq[j+_start]) )
            {
                _matrix[i + j*_N] = 1;
            }
            else
            {
                _matrix[i + j*_N] = 0;
            }
        }
 } 

 //here we update the location of the start and end of the sequence

 if (this->_stop ==  this->_seq.length() - 1)
 {
     return -1;
 }    
 else
 {
    return shift;
 }
 
}
//--------------------------------------------------------------------------------------------------------------------------------------
int AssignSequences::get_maxlen_from_matrix(int &sA_start, int& sA_end, int& sB_start, int& sB_end)
{
    int maxlen = 0;
    int res_i = 0;
    int res_j = 0;
    for(int i = 0; i < _N; i++)
    {
        for(int j = i; j < _N; j++)
        {
            if (_matrix[i + j*_N])
            {
                int sublen = 1;
                for(int k = 1; i-k >= 0 && k+j < _N; k++)
                {
                    if( _matrix[ (i-k) + (j+k)*_N])
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
    
    sA_start  = res_i-maxlen+1 + this->_start;
    sA_end    = res_i + this->_start;
    sB_start  = res_j + this->_start;
    sB_end    = res_j + maxlen - 1 + this->_start;
    //delete [] matrix;
    return maxlen;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bool AssignSequences::are_compatible(char a, char b)
{
    a = toupper(a);
    b = toupper(b);
    if (a == 'U')
    {
        a = 'T';
    }
    if (b == 'U')
    {
        b = 'T';
    }
    if ( (a == 'C' && b == 'G')  ||  (a == 'G' && b == 'C')  )
    {
        return true;
    }
    if( (a == 'A' && b == 'T') ||  (a == 'T' && b == 'A'))
    {
        return true;
    }
    if(this->_wobble &&   ((a == 'G' && b == 'T') ||  (a == 'T' && b == 'G')) )
    {
        return true;
    }

    return false;
}
//--------------------------------------------------------------------------------------------------------------------------------------
double AssignSequences::get_xds(const string& myseq, int start, int end, int maxlen,double c)
{
    //double alpha = 0;

    int l = end-start + 1; //myseq.length();
    double cs  = 0;
    double gs = 0;  //std::count(myseq.begin(), myseq.end(), 'G' ) / double(l);
    double us = 0;  //(std::count(myseq.begin(), myseq.end(), 'T' )  + std::count(myseq.begin(), myseq.end(), 'U' )  )/ double(l);
    double aas = 0; // std::count(myseq.begin(), myseq.end(), 'A' )  / double(l);
    double unknown = 0;
    for (int i = start; i < end; i++)
    {
        switch (toupper(myseq[i]))
        {
            case 'C':  cs += 1; break;
            case 'G':  gs += 1; break;
            case 'T':  us += 1; break;
            case 'U':  us += 1; break;
            case 'A':  aas += 1; break;
            default: unknown += 1; break;
        }
    }
    
    cs   /=  double(l);
    gs   /=  double(l);
    us   /=  double(l);
    aas  /=  double(l);
    unknown /= double(l);

    double alpha = cs*gs + aas*us + gs*cs + gs*us + us*gs + us*aas;

    //alpha = get_alpha(seq)
    double L = double(l);
    double K = maxlen;
    
    double myx =  (-2. * log(L) / (K - c)  ) - log(alpha);

    //#myx = -np.log(alpha) - np.log(L/(2.*a))/maxl
    if(unknown >= 0.6)
    {
      return  std::numeric_limits<double>::quiet_NaN();   
    }
    return myx;
    
}
//--------------------------------------------------------------------------------------------------------------------------------------

int AssignSequences::find_longest_stretch(const string& seq, int &sA_start, int& sA_end, int& sB_start, int& sB_end, bool include_wobble)
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
          if (are_compatible(seq[i],seq[j]))
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
//--------------------------------------------------------------------------------------------------------------------------------------
double AssignSequences::extract_xds_from_matrix(int &maxlen,int& SaS, int& SaE, int& SbS, int& SbE)
{
     //int SaS,SaE,SbS, SbE;

     maxlen = this->get_maxlen_from_matrix(SaS,SaE,SbS, SbE);
     double xds = AssignSequences::get_xds(this->_seq,this->_start,this->_stop,maxlen);
     return xds;
}
//--------------------------------------------------------------------------------------------------------------------------------------
void AssignSequences::scan_and_do_xds(ostream &out)
{
    /*
  this->fill_matrix(0,this->_window_size);
  int total_l = this->_seq.length();
  double xds = extract_xds_from_matrix();
  out << 0 << " " <<  xds << endl;
  for(int i = 0+this->_shift_size; i < total_l; i += this->_shift_size)
  {
      this->shift_and_fill_matrix(i,this->_shift_size);
      xds = extract_xds_from_matrix();
      out << i << " " <<  xds << endl;
  }
  */
  this->fill_matrix(0,_window_size);
  int total_l = _seq.length();
  int SaS,SaE,SbS, SbE;

  int maxlen = 0;// = this->get_maxlen_from_matrix(SaS,SaE,SbS, SbE);
  double xds = extract_xds_from_matrix(maxlen,SaS,SaE,SbS, SbE);
  
  //out << 0 << " " << 0 + _window_size << " " <<  maxlen << " " << xds << " | " << SaS << "  " << SaE << "  " << SbS << " " << SbE  << endl;
  
  out << _contig << "\t" << _cmdline_start+1 << "\t" << _cmdline_start + _window_size << "\t" <<  maxlen << "\t" << xds << "\t|\t " << SaS << "\t" << SaE << "\t" << SbS << "\t" << SbE  << endl;

  for(int i = 0+_shift_size; i <= total_l  - _window_size; i += _shift_size)
  {
      this->shift_and_fill_matrix(i,_shift_size);
      xds = extract_xds_from_matrix(maxlen,SaS,SaE,SbS, SbE);
  
      out << _contig << "\t" << i+_cmdline_start+1 << "\t" << i + _cmdline_start + _window_size << "\t" <<  maxlen << "\t" << xds << "\t|\t " << SaS << "\t" << SaE << "\t" << SbS << "\t" << SbE  << endl;
     
  }


}
//-------------------------------------------------------------------------------------------------------------------------------------
bool AssignSequences::check_segments(string &A, string &B)
{
    if (A.length() != B.length())
    {
        return false;
    }
    for (int i = 0; i < A.length(); i++)
    {
        if ( ! are_compatible(A[i],B[B.length()-1-i]) )
        {
            return false;
        }
    }
    return true;
}

//--------------------------------------------------------------------------------------------------------------------------------------
void AssignSequences::scan_and_do_maxlen(ostream &out)
{
  this->fill_matrix(0,this->_window_size);
  //cout << "Matrix filled " << endl;
  int total_l = this->_seq.length();
  //cout << "Sequence length " << total_l << endl;
  //double xds = extract_xds_from_matrix();
  int SaS,SaE,SbS, SbE;

  int maxlen = this->get_maxlen_from_matrix(SaS,SaE,SbS, SbE);
  
  //out << 0 << ":: " <<  maxlen << "| " << SaS << "  " << SaE << "  " << SbS << " " << SbE  << endl;
  out << 0 << " " <<  maxlen << "| " << SaS << "  " << SaE << "  " << SbS << " " << SbE  << endl;

  for(int i = 0+this->_shift_size; i <= total_l  - _window_size; i += this->_shift_size)
  {
      //cout << " Calling shift and fill " << this->_shift_size << endl;
      this->shift_and_fill_matrix(i,this->_shift_size);
      //this->fill_matrix(i,i+this->_window_size);
      //this->_start += _shift_size;
      int maxlen = this->get_maxlen_from_matrix(SaS,SaE,SbS, SbE);
      //out << i << ":: " <<  maxlen << ": " << SaS << "  " << SaE << "  " << SbS << " " << SbE  << endl;
     
     
      out << i << " " << i+_window_size << " | " <<  maxlen <<  ": " << SaS << "  " << SaE << "  " << SbS << " " << SbE  << " " << endl;
      //out << _seq.substr(SaS,maxlen) << " "  << _seq.substr(SbS,maxlen) << endl; 
  
      //xds = extract_xds_from_matrix();
      //out << i << " " <<  xds << endl;
  }
}

//--------------------------------------------------------------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------------------------------------------------------------
char reverse(char x)
{
    switch(toupper(x))
    {
        case 'A': return 'T';
        case 'T' : return 'A';
        case 'C' : return 'G';
        case 'G' : return 'C';
        default: return 'N'; 
    }
}
//--------------------------------------------------------------------------------------------------------------------------------------
string reverse_complement(const string &s,bool revert=false)
{
    string result = s;
    for(int i = 0; i < s.length(); i++)
    {
        if (revert)
        {
            result[s.length()-1-i] = reverse(s[i]);
        }
        else
        {
            result[i] = reverse(s[i]);
        }
        
    }
    return result;
}
//------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
void process_genome(const string& chrom_seq,const string& output_file_name,const string& chromosome_name,int start,int window_size,int slide_size, string& chrom)
{
    
    //string substr;
    int maxlen = 0;
    int maxlen_r = 0;
    int sAs,sAe,sBs,sBe;
    ofstream out_pos(output_file_name+".positive_sense.xds");
    //ofstream out_neg(output_file_name+".negative_sense");
    //string new_seq = chrom_seq.substr(start,stop-start+1)
    AssignSequences as(chrom_seq,window_size,slide_size,true,start, chrom);
    as.scan_and_do_xds(out_pos);
    
    /*
    int report = 0; 
    for (int i = start; i < stop; i += slide_size)
    {
        //maxlen = 0;
        //cout << "Processing " << i << endl;
        string substr = chrom_seq.substr(i,window_size);
        string revsub = reverse_complement(substr);
        //cout << "Substr " << substr.length() << endl;
        maxlen = find_longest_stretch(substr,sAs,sAe,sBs,sBe);
        if (maxlen >= length_cut)
        {
            out_pos << maxlen << "\t" << chromosome_name <<  " " << i+sAs << " " << i+sAe << " " << i+sBs << " " << i+sBe << endl;
        }
        maxlen_r = find_longest_stretch(revsub,sAs,sAe,sBs,sBe);
        if (maxlen_r >= length_cut)
        {
            out_neg << maxlen_r << "\t" << chromosome_name <<  " " << i+sAs << " " << i+sAe << " " << i+sBs << " " << i+sBe << endl;
        }

        if( int(double(100.*i)/stop) >= report+10 )
        {
            report = int(double(100.*i)/stop);
            cerr << "Processed " <<  int(double(100.*i)/stop) << "% of chromosome " << chromosome_name << endl;
        } 
    }
    */

}

//--------------------------------------------------------------------------------------------------------------------------------------
string load_chromosome(const string hg38_file,const string& chromosome)
{
    std::ifstream input(hg38_file);
    
    list<name_seq> loaded_seqs;
    if(!input.good()){
        std::cerr << "Error opening '"<< hg38_file <<"'. Bailing out." << std::endl;
        throw runtime_error("Could not open file");
    }
 
    std::string line, name, content;
    while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Print out what we read from the last entry
                //std::cout << name << ": " << content << std::endl;
                //std::cout << " looking for:" << chromosome << ";" << endl;
                //loaded_seqs.push_back( make_pair(name,content)); 
                if(trim(name) == chromosome) 
                { 
                    //std::cout << " Loaded " << chromosome << ";" << endl;
                    transform(content.begin(),content.end(), content.begin(),::toupper);
                     return content;

                }
                name.clear();


            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }
    if( !name.empty() ){ // Print out what we read from the last entry
        //std::cout << name << " : " << content << std::endl;
        transform(content.begin(),content.end(), content.begin(),::toupper);
        //loaded_seqs.push_back( make_pair(name,content));
        //cout << "loaded " << name << "|" << endl;
        if ( trim(name) == chromosome )
            return content;
    }

    //cout << "Loaded " << loaded_seqs.size() << " sequences " << endl;
    //loaded_seqs.push_back( make_pair(name,content));
    throw runtime_error(string("Could not find chromosome") + chromosome);
    return content;

}
//---------------------------------------------------------------------------------------------------------------------------------------

/* void process_sequence_file(const char *fastafile,int window_size,int step_size,string outputfilename)
{
       list<name_seq> rs = read_fasta(string(fastafile));
       cerr << "Loaded " << rs.size() << " sequences " << endl;

       
       for (auto i = rs.begin(); i != rs.end(); ++i)
       {
           //cerr << "Processing: Found: " << (*i).first << " - -- "  << (*i).second << endl;
           //cout << (*i).first << " "  << endl;
           //ProcessSeq ps((*i).second,6,(*i).second.length()/2,(*i).first);
           //cout << "Object constructed" << endl;
           int sAs,sAe,sBs,sBe;
           //int maxlen = find_longest_stretch((*i).second,sAs,sAe,sBs,sBe);
           AssignSequences as((*i).second, window_size,step_size);
           //as.scan_and_do_maxlen(cout);
           as.scan_and_do_xds(cout);
           //cout << (*i).first << " " <<  maxlen <<  " " << sAs << "  "  << sAe  << " " << sBs << "  "  << sBe <<  endl;
           //cout << (*i).second.length() << " "  << maxlen << " " <<  AssignSequences::get_xds((*i).second,maxlen) << endl;
       }
 
} */

//---------------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
    
    if(argc != 8)
	{
		cerr << "Usage: " << argv[0] << " fasta_file_with_genome chromosome_name start-1-based stop shiftsize sequence_length outputfilename " << endl;
		cerr << "Output: contig window_start window_end length_of_complem_segment double_stranded_force | start_of_left_segment end_of_right_segment start_of_right_segment end_of_right_segment" << endl;
		 cerr << "note: if start or stop parameters is smaller than 0, the program scans the entire chromosome" << endl; 
        return 1;
	}
    try
    {
       //cout << reverse_complement(string("AAAATTC")) << endl;
   
       string genome_file(argv[1]);
       string chrom(argv[2]);
       int start = atoi(argv[3]) - 1; // convert 1-based to 0-based
       int stop  = atoi(argv[4]);  // no need to convert since 0-based
           // assumes half-open interval and 1-based assumes open interval
       //int minlength = atoi(argv[4]);
       int shiftsize = atoi(argv[5]);
       int seqlen    = atoi(argv[6]);
       string outputfname(argv[7]);
       string chromseq = load_chromosome(genome_file,chrom);
       if (start < 0 || stop < 0)
       {
           start = 0;
           stop = chromseq.length();
       }
       else
       {
           chromseq = chromseq.substr(start,stop-start+1);
       }
       cerr << "Loaded chromosome seq " << chromseq.length() << " long " << endl;
       //cerr << chromseq << endl;
       process_genome(chromseq,outputfname,chrom,start,seqlen,shiftsize,chrom);


    }
    catch(exception& e)
    {
        cerr << e.what() << endl;
        cerr << " finished with error" << endl;
        return -1;
    }
    return 0;


}
//--------------------------------------------------------------------------------------------------------------------------------------
