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


bool are_compatible(char a, char b,bool include_wobble=true)
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
    if(include_wobble &&   ((a == 'G' && b == 'T') ||  (a == 'T' && b == 'G')) )
    {
        return true;
    }

    return false;
}
//--------------------------------------------------------------------------------------------------------------------------------------
double get_xds(const string& myseq, int maxlen,double c = -2.2)
{
    //double alpha = 0;

    int l = myseq.length();
    double cs = std::count(myseq.begin(), myseq.end(), 'C' ) / double(l);
    double gs = std::count(myseq.begin(), myseq.end(), 'G' ) / double(l);
    double us = (std::count(myseq.begin(), myseq.end(), 'T' )  + std::count(myseq.begin(), myseq.end(), 'U' )  )/ double(l);
    double aas = std::count(myseq.begin(), myseq.end(), 'A' )  / double(l);
    double alpha = cs*gs + aas*us + gs*cs + gs*us + us*gs + us*aas;

    //alpha = get_alpha(seq)
    double L = double(l);
    double K = maxlen;
    
    double myx =  (-2. * log(L) / (K - c)  ) - log(alpha);

    //#myx = -np.log(alpha) - np.log(L/(2.*a))/maxl
    return myx;
}
//--------------------------------------------------------------------------------------------------------------------------------------

int find_longest_stretch(const string& seq, int &sA_start, int& sA_end, int& sB_start, int& sB_end, bool include_wobble=true)
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
void process_genome(const string& chrom_seq,const string& output_file_name,const string& chromosome_name,int start,int stop,int window_size,int slide_size,int length_cut)
{
    //string substr;
    int maxlen = 0;
    int maxlen_r = 0;
    int sAs,sAe,sBs,sBe;
    ofstream out_pos(output_file_name+".positive_sense");
    ofstream out_neg(output_file_name+".negative_sense");
    
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
                //std::cout << name << " : " << content << std::endl;
                //loaded_seqs.push_back( make_pair(name,content)); 
                if(name == chromosome) 
                { 
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
        if ( name == chromosome )
            return content;
    }

    //cout << "Loaded " << loaded_seqs.size() << " sequences " << endl;
    //loaded_seqs.push_back( make_pair(name,content));
    throw runtime_error(string("Could not find chromosome") + chromosome);
    return content;

}
//---------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------
int genome_main(int argc, char **argv)
{
    string genome_file = "/media/petr/SCRATCH1/vochomurka/home_dir/petr/projects/CpG_evolution/hg38.fa";
    //string genome_file = "/tmp/example.fa";
    
    if(argc != 8)
	{
		cerr << "Usage: " << argv[0] << " chromosome start stop minlength shiftsize sequence_length outputfilename " << endl;
		return 1;
	}
    try
    {
       //cout << reverse_complement(string("AAAATTC")) << endl;
       string chrom(argv[1]);
       int start = atoi(argv[2]);
       int stop  = atoi(argv[3]);
       int minlength = atoi(argv[4]);
       int shiftsize = atoi(argv[5]);
       int seqlen    = atoi(argv[6]);
       string outputfname(argv[7]);
       string chromseq = load_chromosome(genome_file,chrom);
       if (start < 0 || stop < 0)
       {
           start = 0;
           stop = chromseq.length();
       }
       cerr << "Loaded chromosome seq " << chromseq.length() << " long " << endl;
       process_genome(chromseq,outputfname,chrom,start,stop,seqlen,shiftsize, minlength);


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
//---------------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
    if(argc != 2)
	{
		cerr << "Usage: " << argv[0] << " fastafile " << endl;

		cerr << "Outputs data in the following format (in the order of the sequences in the fasta file): sequence_length max_complementary_segment_length double_stranded_force "  << endl;
		return 1;
	}
    try
    {
       list<name_seq> rs = read_fasta(argv[1]);
       cerr << "Loaded " << rs.size() << " sequences " << endl;

       for (auto i = rs.begin(); i != rs.end(); ++i)
       {
           //cerr << "Processing: Found: " << (*i).first << " - -- "  << (*i).second << endl;
           //cout << (*i).first << " "  << endl;
           //ProcessSeq ps((*i).second,6,(*i).second.length()/2,(*i).first);
           //cout << "Object constructed" << endl;
           int sAs,sAe,sBs,sBe;
           int maxlen = find_longest_stretch((*i).second,sAs,sAe,sBs,sBe);
           //cout << (*i).first << " " <<  maxlen <<  " " << sAs << "  "  << sAe  << " " << sBs << "  "  << sBe <<  endl;
           cout << (*i).second.length() << " "  << maxlen << " " <<  get_xds((*i).second,maxlen) << endl;
       }
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
