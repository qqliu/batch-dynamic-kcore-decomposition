
# include <iostream>
# include <set>
# include <map>
# include <utility>
# include <fstream>
# include <string>
# include <sstream>

using namespace std;

void getUndirectGraph(void);

int main(void)
{
	getUndirectGraph();



	return 0;
}

void getUndirectGraph(void)
{
	set<pair<string, string> > setEdges;
	ifstream fin("./uk-2002.txt");
	if (fin)
	{
		int count = 1;
		string line;
		while (getline(fin, line))
		{
			if (count == 1)
			{
				size_t pos = line.find("\t");
				size_t size = line.size();
				if (pos != string::npos)
				{
					string strNode1 = line.substr(0, pos);
					string strNode2 = line.substr(pos+1, size);

					stringstream ss1;
					int intNode1;
					ss1 << strNode1;
					ss1 >> intNode1;
					stringstream ss2;
					int intNode2;
					ss2 << strNode2;
					ss2 >> intNode2;

					if (intNode1 < intNode2)
					{
						pair<string, string> pairEdge(strNode1, strNode2);
						setEdges.insert(pairEdge);
					}
					else if (intNode1 > intNode2)
					{
						pair<string, string> pairEdge(strNode2, strNode1);
						setEdges.insert(pairEdge);
					}
					else
					{
						;
					}
				}
			}
			else
			{
				count += 1;
			}
		}
	}
	else
	{
		cout << "input error" << endl;
	}
	fin.close();

	ofstream fout("./uk_undirect.txt");
	if (fout)
	{
		for (set<pair<string, string> >::iterator it = setEdges.begin(); it != setEdges.end(); it++)
		{
			//cout << it->first << "\t" << it->second << endl;
			string strContent = it->first + "\t" + it->second;
			fout << strContent << endl;
		}
	}
	else
	{
		cout << "output error" << endl;
	}
	fout.close();
}
