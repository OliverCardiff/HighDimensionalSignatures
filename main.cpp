#include <stdio.h>
#include <getopt.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "Node.hpp"
#include "WedgeMatrix.hpp"
#include "TreeManager.hpp"
#include "ThreadManager.hpp"

using namespace std;

int ThreadManager::SUB_THRESHOLD = 500;

ThreadManager * Tree::TMAN = 0;

float WedgeMatrix::PLOIDY = 2;
unsigned int WedgeMatrix::SETSIZE = 4;

unsigned int Node::PLOIDY = 2;
unsigned int Node::DEPTH = 25;
unsigned short Node::NMAX = 2;
unsigned int Node::SETSIZE = 4;
unsigned int Node::AGGS = 0;
vector<char> * Node::CHARSET;
char * Node::PREFIX = 0;

//unsigned int NetNode::DEPTH = 25;
//unsigned int NetNode::SETSIZE = 4;
//vector<char> * NetNode::CHARSET;

vector<unsigned int> Node::TERMINUS;

int main(int argc, char **argv)
{
	char *genome = 0;
	char *prefix = 0;
	char *depthFile = 0;
	int depth = 25;
	int ploidy = 2;
	int c;
	int threads = 4;
	short nmax = 1;
	
	bool hasGenome = false;
	bool doSubs = true;
	bool hasDepths = false;
	bool hasPrefix = false; 
	bool wants_help = false;
	bool doNull = true;
	bool goPro = true;
	
	while((c = getopt (argc, argv, "x:n:t:g:d:o:p:hsu")) != -1)
	{
		switch(c)
		{
			case 'g':
				hasGenome = true;
				genome = optarg;
				break;
			case 's':
				doSubs = false;
				break;
			case 'u':
				doNull = false;
				break;
			case 'x':
				hasDepths = true;
				depthFile = optarg;
				break;
			case 't':
				threads = atoi(optarg);
				break;
			case 'n':
				nmax = atoi(optarg);
				break;
			case 'o':
				hasPrefix = true;
				prefix = optarg;
				break;
			case 'h':
				wants_help = true;
				break;
			case 'd':
				depth = atoi(optarg);
				break;
			case 'p':
				ploidy = atoi(optarg);;
				break;
			case '?':
				if(optopt == '1' || optopt == '2' || optopt == 'o'|| optopt == 'c')
				{
					fprintf(stderr, "Option -%c requires an argument.\n", optopt);
				}
				else if(isprint (optopt))
				{
					fprintf(stderr, "Unknown Option '-%c'.\n", optopt);
				}
				else
				{
					fprintf(stderr, "Unknown Option '\\x%x'.\n", optopt);
				}
				return 1;
			default:
				abort();
		}
	}
	
	if(wants_help)
	{
		cout << "\nUGPep: Protein Structural Complexity Assessment Tool\n\n";
		cout << "Required Input:\n";
		cout << "-g\t\t\tThe fasta formatted genome file\n";
		cout << "-o\t\t\tThe output file prefix, includes pathway\n\n";
		cout << "Optional Inputs:\n";
		cout << "-p\t\t\tThe genomic ploidy of the organism <int> [2]\n";
		cout << "-x\t\t\tPaired read-depth file\n";
		cout << "-d\t\t\tMax tree depth/k-mer length <int> [25]\n";
		cout << "-t\t\t\tMax threads to run (2*cpus reccommended) <int>[4]\n";
		cout << "-n\t\t\tNumber of Ns in kmers <int>[2]\n";
		cout << "-s\t\t\tFlag, prevents escaped frequency subtraction\n";
		cout << "-u\t\t\tFlag, prevents null tree generation\n";
		cout << "-h\t\t\tPrints this screen\n\n";
		cout << "Typical Usage:\n\n";
		cout << "./UGPep -g my_proteome.pep -d 28 -o ~/file/path/to/my_kd_out\n\n";
		return 0;
	}
	
	if(!goPro)
	{
		char nchr[] = {'A', 'T', 'G', 'C'};
		Node::SETSIZE = 4;
		//NetNode::SETSIZE = 4;
		WedgeMatrix::SETSIZE = 4;
		Node::NMAX = nmax;
		Node::CHARSET = new vector<char>(nchr, nchr + sizeof(nchr)/sizeof(char));
		//NetNode::CHARSET = new vector<char>(nchr, nchr + sizeof(nchr)/sizeof(char));
	}
	else
	{
		char nchr[] = {'A', 'R', 'N', 'D', 'X', 'C', 'E', 'Q', 'Z', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
		Node::SETSIZE = 22;
		//NetNode::SETSIZE = 22;
		WedgeMatrix::SETSIZE = 22;
		Node::NMAX = nmax;
		Node::CHARSET = new vector<char>(nchr, nchr + sizeof(nchr)/sizeof(char));
		//NetNode::CHARSET = new vector<char>(nchr, nchr + sizeof(nchr)/sizeof(char));
	}
	
	if(!hasGenome)
	{
		cout << "You need to specify the genome fasta file with '-g' <genome.fa>\n";
	}
	if(!hasPrefix)
	{
		cout << "You need to specify the output prefix with '-o' <prefix>\n";
	}
	
	if(hasGenome && hasPrefix)
	{
		Node::PLOIDY = ploidy;
		WedgeMatrix::PLOIDY = (float)ploidy;
		Node::DEPTH = depth;
		Node::PREFIX = prefix;
		Node::TERMINUS.resize(10000);
		
		Tree::TMAN = new ThreadManager(threads);
		TreeManager tr(depth, doSubs, doNull);
		
		/*Node boy(0,false);
		map<unsigned int, unsigned int> mp;
		Node * pboy = 0;
		vector<Node*> vec;
		unsigned int tint = 0;
		pair<unsigned int, unsigned int> ipair(0,0);
		
		cout << "Nsize: " << sizeof(boy) << endl;
		cout << "Msize: " << sizeof(mp) << endl;
		cout << "Psize: " << sizeof(pboy) << endl;
		cout << "uint : " << sizeof(tint) << endl;
		cout << "vsize: " << sizeof(vec) << endl;
		cout << "Ipair: " << sizeof(ipair) << endl;*/
		
		
		tr.RunZDSCS(genome, hasDepths, depthFile);
		//tr.ZWedgeSummary();
		tr.WriteWedgeFiles(prefix, false);
		tr.SubtractAndResolve();
		if(doNull)
		{
			tr.WriteWedgeFiles(prefix, true);
		}
		tr.CollapseSDs();
		tr.WriteSDMatrices(prefix);
		//string boi;
		//cin >> boi; 
	}
	
	return 0;
}
