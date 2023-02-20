// Copyright (C) 2003, 2004 by BiRC -- Bioinformatics Research Center
//                             University of Aarhus, Denmark
//                             Contact: Thomas Mailund <mailund@birc.dk>

#include "label-map.hh"
using namespace std;

size_t
LabelMap::push(string label) throw(AlreadyPushedEx)
{
    if (_map.find(label) != _map.end()) throw AlreadyPushedEx(label);
    _names.push_back(label);
    return _map[label] = _count++;
}

size_t
LabelMap::operator[](string label) const throw(UnkownLabelEx)
{
    map<string, size_t>::const_iterator i = _map.find(label);
    if (i == _map.end()) throw UnkownLabelEx(label);
    return i->second;
}

std::string
LabelMap::name(unsigned int idx) const throw(std::out_of_range)
{
    return _names.at(idx);
}
