#include "Tree.hpp"
#include <fstream>
#include "Node.hpp"
#include <iostream>
#include <math.h>
#include <ctime>

using namespace std;

Tree::Tree(int depth, bool IsNull, bool doSubs)
{
	_hasDepths = false;
	_subtract = doSubs;
	_depth = depth;
	_null = IsNull;
	_root = new Node(0, true);
	_zWedge = new WedgeMatrix(depth, true);
	_ordinance = 0;
	_skipRate = 0;
	_maxPenalty = 0;
	_psRate = 0;
	srand(time(0));
}

Tree::~Tree()
{
	//delete _zWedge;
	delete _root;
}

void Tree::SetSkip(double skip)
{
	_skipRate = skip;
}

WedgeMatrix * Tree::GetZWedge()
{
	return _zWedge;
}

double Tree::GetSEscape()
{
	return _psRate;
}

double Tree::GetHeadFreq()
{
	return _root->GetFreq();
}

void Tree::ResolveParetos()
{
	_zWedge->ResolveParetos();
}

void Tree::SetNullReference(WedgeMatrix * nullMat)
{
	_zWedge->SetNull(nullMat);
}

void Tree::ExpandTree()
{
	vector<bool> npath;
	
	for(int i = 0; i < 10000; i++)
	{
		Node::TERMINUS[i] = 0;
	}
	
	if(_subtract)
	{
		_maxPenalty = _root->BackPropSubtract(0);
	}

	_zWedge->AddHeadFreq(_root->GetFreq());
	TMAN->SetZWedge(_zWedge);
	
	_root->FindNextExpand(0, &npath, 0,0,0, _zWedge, 0);
	
	_ordinance = Node::UnexplodedOrdinance();
	if(_subtract)
	{
		_psRate = Node::PairedSeqsRate();
	}
	else
	{
		_psRate = 0;
	}
	_plRate = Node::PloidyRate();
	_escapeS = Node::TerminalStruct();
} 

void Tree::FirstTreeSummary()
{
	for(int i = 0; i < 10000; i++)
	{
		Node::TERMINUS[i] = 0;
	}
	
	unsigned int rFreq = _root->GetFreq();
	
	if(_subtract)
	{
		_maxPenalty = _root->BackPropSubtract(0);
	}
	
	cout << "Prior Input Frequency: \t" << rFreq << "\n";
	cout << "Escaped Freq Penalty: \t" << _maxPenalty << "\n";
	cout << "Tree thinned by: \t" << ((double)_maxPenalty / (double)rFreq ) * 100 << " %" << endl;
	
	_ordinance = Node::UnexplodedOrdinance();
	if(_subtract)
	{
		_psRate = Node::PairedSeqsRate();
	}
	else
	{
		_psRate = 0;
	}
	_plRate = Node::PloidyRate();
	_escapeS = Node::TerminalStruct();
}

void Tree::StateZSummary(bool escapes)
{
	cout << "Normalised S =\t" << _zWedge->GetNormalisedS() << "\n";
	cout << "Set Pareto =  \t" << _zWedge->GetSumP() << "\n";
	cout << "Pareto N =  \t" << _zWedge->GetParetoN() << "\n"; 
	
	cout << "Root freq =\t" << _root->GetFreq() << "\n";
	cout << "Actual S  =\t" << _zWedge->GetSumStruct() << "\n";
	cout << "Actual F  =\t" << _zWedge->GetSumFs() << "\n";
	
	double sumS = _zWedge->GetSumStruct();
	
	double bounds = log(_root->GetFreq())/log(Node::SETSIZE);
	
	if(escapes)
	{
		cout << "\nZero->K Entropy: \t\t" << ((1 - _ordinance) * 100) << " %\n";
		cout << "Unused Frequency Rate: \t\t" << (_psRate * 100) << " %\n";
		
		if(Node::PLOIDY > 1)
		{
			cout << "..of which est. by ploidy: \t" << _plRate << " %\n";
			cout << "Min structure escape: \t\t" << (_escapeS/(sumS + _escapeS)) * 100 << " %\n";
		}
	}
	cout << "\nDepth: " << Node::DEPTH << "\n";
	cout << "Alphabet of " << Node::SETSIZE << ", null saturation at power: " << bounds << endl;
	
}

void Tree::ReadDepths(char * depthFile)
{
	string line;
	
	ifstream mygen (depthFile);
	if(mygen.is_open())
	{
		while(getline(mygen, line))
		{
			_parallelDepths.push_back(atoi(line.c_str()));
		}
	}
	mygen.close();
	
	_hasDepths = true;
}

void Tree::ReadFasta(char * genome)
{
	string line;
	string seqs;
	int inputctr = 0;
	int lenctr = 0;
	int cov = 1;
	int ind = 0;
	
	ifstream mygen (genome);
	if(mygen.is_open())
	{
		while(getline(mygen, line))
		{
			const char * cchar = line.c_str();
			if(cchar[0] != '>')
			{
				seqs = seqs + line; 
			}
			else
			{
				if((int)seqs.length() > _depth)
				{
					if(_hasDepths)
					{
						cov = (int)floor(log2((double)_parallelDepths[ind]));
						ind++;
					}
					if(cov > 0)
					{
						if(_null)
						{
							ReadNullString(&seqs, cov);
						}
						else
						{
							ReadInString(&seqs, cov);
						}
						inputctr++;
						lenctr += seqs.length();
					}
				}
				seqs = "";
			}
		}
		mygen.close();
	}
	if((int)seqs.length() > _depth)
	{
		if(_hasDepths)
		{
			cov = (int)floor(log2((double)_parallelDepths[ind]));
			ind++;
		}
		if(cov > 0)
		{
			if(_null)
			{
				ReadNullString(&seqs, cov);
			}
			else
			{
				ReadInString(&seqs, cov);
			}
			inputctr++;
			lenctr += seqs.length();
		}
	}
}

void Tree::GetZeroS()
{
	_root->RunZDFS(_zWedge);
}

void Tree::ReadNullString(string * sts, int cov)
{
	int len = sts->length() - _depth;
	double testSkip = 1;
	
	for(int i = 0; i < len; i++)
	{
		string n2 = sts->substr(i, _depth);
		
		testSkip = (double)rand()/(double)RAND_MAX;
		
		if(testSkip > _skipRate)
		{
			random_shuffle(n2.begin(), n2.end());

			_root->SendString(n2.c_str(), 0, cov);
			if(Node::SETSIZE == 4)
			{
				RevComp(&n2);
				_root->SendString(n2.c_str(), 0, cov);
			}
			else
			{
				reverse(n2.begin(), n2.end());
				_root->SendString(n2.c_str(), 0, cov);
			}
		}
	}
}

void Tree::ReadInString(string * sts, int cov)
{
	int len = sts->length() - _depth;

	for(int i = 0; i < len; i++)
	{
		string n2 = sts->substr(i, _depth);
		
		_root->SendString(n2.c_str(), 0, cov);
		if(Node::SETSIZE == 4)
		{
			RevComp(&n2);
			_root->SendString(n2.c_str(), 0, cov);
		}
		else
		{
			//reverse(n2.begin(), n2.end());
			//_root->SendString(n2.c_str(), 0, cov);
		}
	}
}

void Tree::RevComp(string * sts)
{
	string * next = sts;
	reverse(next->begin(), next->end());
	
	for(unsigned int i = 0; i < next->length(); i++)
	{
		char t = (*next)[i];
		if(t == 'A') (*next)[i] = 'T';
		else if(t == 'T') (*next)[i] = 'A';
		else if(t == 'G') (*next)[i] = 'C';
		else if(t == 'C') (*next)[i] = 'G';
	}
}