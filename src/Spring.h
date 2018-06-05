#pragma once
#ifndef __Spring__
#define __Spring__

#include <memory>

class Node;

class Spring
{
public:
	Spring(std::shared_ptr<Node> p0, std::shared_ptr<Node> p1);
	virtual ~Spring();
	
	std::shared_ptr<Node> p0;
	std::shared_ptr<Node> p1;
	double E;
	double L;
};

#endif
