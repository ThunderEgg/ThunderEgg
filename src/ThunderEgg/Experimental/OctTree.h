/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  ThunderEgg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

#ifndef THUNDEREGG_EXPERIMENTAL_OCTREE_H
#define THUNDEREGG_EXPERIMENTAL_OCTREE_H

#include <ThunderEgg/Experimental/OctNode.h>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>

namespace ThunderEgg
{
namespace Experimental
{
/**
 * @brief Represents an oct-tree.
 */
template <int D> struct Tree {
	/**
	 * @brief Map of node id to node objects.
	 */
	std::map<int, Node<D>> nodes;
	/**
	 * @brief Map of level to a node object on that level. With 0 begin the root of the tree.
	 */
	std::map<int, Node<D> *> levels;
	/**
	 * @brief The id of the root.
	 */
	int root;
	/**
	 * @brief The number of levels in this tree.
	 */
	int num_levels;
	/**
	 * @brief The maximum id value used in this tree.
	 */
	int max_id;
	/**
	 * @brief Create new tree with only a root node.
	 */
	Tree();
	/**
	 * @brief Read an Tree from a file.
	 *
	 * @param file_name the file to read from.
	 */
	Tree(std::string file_name);
	/**
	 * @brief Refine the tree to add one more level.
	 */
	void refineLeaves();
	/**
	 * @brief Refine a particular node.
	 *
	 * @param n the node to refine.
	 */
	void refineNode(Node<D> &n);
	void zoltanBalanceDomains();
};
template <int D> inline Tree<D>::Tree()
{
	Node<D> new_root;
	new_root.id = 0;
	new_root.lengths.fill(1);
	new_root.starts.fill(0);
	new_root.level = 0;
	nodes[0]       = new_root;
	levels[1]      = &nodes[0];
	root           = 0;
	max_id         = 0;
	num_levels     = 1;
}
template <int D> inline Tree<D>::Tree(std::string file_name)
{
	using namespace std;
	ifstream input(file_name, ios_base::binary);
	int      num_nodes;
	int      num_trees;
	input.read((char *) &num_nodes, 4);
	input.read((char *) &num_trees, 4);
	num_levels = 0;
	max_id     = 0;
	for (int i = 0; i < num_nodes; i++) {
		Node<D> n;
		// id level parent
		input.read((char *) &n.id, 3 * 4);
		// length and starts
		input.read((char *) &n.lengths[0], sizeof(n.lengths));
		input.read((char *) &n.starts[0], sizeof(n.starts));

		input.read((char *) &n.nbr_id[0], sizeof(n.nbr_id));
		input.read((char *) &n.child_id[0], sizeof(n.child_id));

		if (i == 0) {
			root = n.id;
		}
		max_id      = max(max_id, n.id);
		nodes[n.id] = n;

		if (n.level > num_levels) {
			num_levels = n.level;
		}
		levels[n.level] = &nodes[n.id];
	}
}
template <int D> inline void Tree<D>::refineLeaves()
{
	using namespace std;
	Node<D> root_node = nodes[root];
	Node<D> child     = root_node;
	int     level     = 0;
	while (child.hasChildren()) {
		child = nodes[child.child_id[0]];
		level++;
	}
	deque<pair<int, int>> q;
	set<pair<int, int>>   qed;
	q.push_back(make_pair(level, child.id));
	qed.insert(make_pair(level, child.id));

	while (!q.empty()) {
		pair<int, int> p     = q.front();
		Node<D>        n     = nodes[p.second];
		int            level = p.first;
		q.pop_front();

		// set and enqueue nbrs
		for (Side<D> s : Side<D>::getValues()) {
			// coarser
			if (n.nbrId(s) == -1 && n.parent != -1 && nodes[n.parent].nbrId(s) != -1) {
				Node<D>        parent = nodes[n.parent];
				Node<D>        nbr    = nodes[parent.nbrId(s)];
				pair<int, int> p(level - 1, nbr.id);
				if (!qed.count(p)) {
					q.push_back(p);
					qed.insert(p);
				}
				// finer
			} else if (n.nbrId(s) != -1 && nodes[n.nbrId(s)].hasChildren()) {
				Node<D> nbr  = nodes[n.nbrId(s)];
				auto    octs = Orthant<D>::getValuesOnSide(s.opposite());
				for (size_t i = 0; i < octs.size(); i++) {
					int            id = nbr.childId(octs[i]);
					pair<int, int> p(level + 1, id);
					if (!qed.count(p)) {
						q.push_back(p);
						qed.insert(p);
					}
				}
				// normal
			} else if (n.nbrId(s) != -1) {
				int            id = n.nbrId(s);
				pair<int, int> p(level, id);
				if (!qed.count(p)) {
					q.push_back(p);
					qed.insert(p);
				}
			}
		}
	}
	for (auto p : qed) {
		refineNode(nodes[p.second]);
	}
	levels[num_levels + 1] = &nodes.at(levels[num_levels]->child_id[0]);
	this->num_levels++;
}
template <int D> inline void Tree<D>::refineNode(Node<D> &n)
{
	std::array<Node<D>, Orthant<D>::num_orthants> new_children;
	for (Orthant<D> o : Orthant<D>::getValues()) {
		new_children[o.getIndex()] = Node<D>(n, o);
		Node<D> &child             = new_children[o.getIndex()];
		max_id++;
		child.id                 = max_id;
		n.child_id[o.getIndex()] = max_id;
	}
	// set new neighbors
	for (Orthant<D> o : Orthant<D>::getValues()) {
		for (Side<D> s : o.getInteriorSides()) {
			new_children[o.getIndex()].nbrId(s) = new_children[o.getNbrOnSide(s).getIndex()].id;
		}
	}

	// set outer neighbors
	for (Side<D> s : Side<D>::getValues()) {
		if (n.hasNbr(s) && nodes[n.nbrId(s)].hasChildren()) {
			Node<D> &nbr = nodes[n.nbrId(s)];
			for (Orthant<D> o : Orthant<D>::getValuesOnSide(s)) {
				Node<D> &child                = new_children[o.getIndex()];
				Node<D> &nbr_child            = nodes[nbr.childId(o.getNbrOnSide(s))];
				child.nbrId(s)                = nbr_child.id;
				nbr_child.nbrId(s.opposite()) = child.id;
			}
		}
	}
	// add nodes
	for (Node<D> child : new_children) {
		nodes[child.id] = child;
	}
}
} // namespace Experimental
} // namespace ThunderEgg
#endif
