#include <cstdio>
#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <utility>
#include <numeric>
#include <unordered_map>
#include <string>
#include <algorithm>

using namespace std;

int main(int argc, char *argv[]) {
	if (argc != 4) {
		fprintf(stderr,"usage: %s <graph-file> <vertex-number> <edge-number>\n", argv[0]);
		return -1;	
	}
	std::set<std::pair<int,int>> edges_set;
	std::set<int> vertices_set;
	std::vector<std::pair<int,int>> edges_list;

	std::fstream fin;
	fin.open(std::string(argv[1]) + ".txt",std::ios_base::in);
	if (fin.fail()) {
		fprintf(stderr, "open %s failed.\n", argv[1]);
		return -1;
	}
	std::unordered_map<int,std::vector<int>> graph;
	int a,b;
	while (fin >> a >> b) {
		if (a >= b) continue;
		edges_set.insert(std::make_pair(a,b));
		auto it = graph.find(a);
		if (it != graph.end()) {
			it->second.push_back(b);
		} else {
			graph[a] = std::vector<int>{b};
		}

		it = graph.find(b);
		if (it != graph.end()) {
			it->second.push_back(a);
		} else {
			graph[b] = std::vector<int>{a};
		}

		vertices_set.insert(a);
		vertices_set.insert(b);
	}

	std::cout << "graph read completed" << std::endl;
	edges_list.insert(edges_list.begin(),edges_set.begin(),edges_set.end());
	edges_set.clear();
	

	auto vertices_size = vertices_set.size();
	std::cout << "vertices size: " << vertices_size << std::endl;
	auto size = edges_list.size();
	std::vector<int> vec,node_vec;
	vec.resize(size);
	std::iota(vec.begin(),vec.end(),0);
	node_vec.resize(vertices_size);
	std::iota(node_vec.begin(),node_vec.end(),0);

	for (auto c = 0; c < 1; c++) {
		std::fstream fout;
		char buf[64];
		sprintf(buf,"%s_syn_%d.txt",argv[1],c);
		fout.open(buf,std::ios_base::in | std::ios_base::out | std::ios_base::trunc);
		if (fout.fail()) {
			fprintf(stderr,"open %s failed.\n", buf);
			return -1;
		}

		std::random_shuffle(node_vec.begin(), node_vec.end());
		auto select_vertices_nums = std::stoi(argv[2]);
		std::set<pair<int,int>> select_edges;

		std::cout << graph.size() << "\t" << node_vec.size() << std::endl;
		// select vertices
		for (auto i = 0; i < select_vertices_nums; ) {
			auto node = node_vec[i];
			auto it = graph.find(node);
			if (it != graph.end()) {
				for (auto j = 0; j < it->second.size(); j++) {
					select_edges.insert(make_pair(node,it->second[j]));
				}
				i++;
			}

		}

		std::cout << "vertices done" << std::endl;

		// select edges
		auto select_edges_nums = std::stoi(argv[3]);
		auto count = select_edges.size();
		auto target = select_edges_nums + count;

		if (target > edges_list.size()) {
			std::cout << "too many edges select!" << std::endl;
			return -1;
		}
		while (count < target) {
			std::random_shuffle(vec.begin(),vec.end());
			
			for (auto i = 0; i < select_edges_nums && count < target; i++) {
				select_edges.insert(edges_list[vec[i]]);
				count = select_edges.size();
			}
		}		

		for (auto e : select_edges) {
			fout << e.first << ' ' << e.second << '\n';
		}

		fout.close();
	}
	
	return 0;
}
