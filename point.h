/*
 * point.h
 *
 *  Created on: Feb 11, 2013
 *      Author: Vadim Fedorov
 */
/**
 * Copyright (C) 2013, Vadim Fedorov <vadim.fedorov@upf.edu>
 *  
 *  This program is free software: you can use, modify and/or
 *  redistribute it under the terms of the simplified BSD
 *  License. You should have received a copy of this license along
 *  this program. If not, see
 *  <http://www.opensource.org/licenses/bsd-license.html>.
 */

#ifndef POINT_H_
#define POINT_H_

struct Point
{
	float x, y;

	Point();
	Point(float x, float y);

	bool operator== (const Point &p) const;
	bool operator!= (const Point &p) const;
	Point& operator= (const Point &p);
	Point& operator+= (const Point &p);
	Point& operator-= (const Point &p);
	const Point operator+ (const Point &p) const;
	const Point operator- (const Point &p) const;

	friend inline bool operator< (const Point& lhs, const Point& rhs);
	friend inline bool operator> (const Point& lhs, const Point& rhs);
	friend inline bool operator<=(const Point& lhs, const Point& rhs);
	friend inline bool operator>=(const Point& lhs, const Point& rhs);

};

// NOTE: definitions are in header in order to overload two argument versions.
inline bool operator< (const Point& lhs, const Point& rhs)
{
	return (lhs.y < rhs.y) || (lhs.y == rhs.y && lhs.x < rhs.x);
}
inline bool operator> (const Point& lhs, const Point& rhs) { return operator< (rhs,lhs); }
inline bool operator<= (const Point& lhs, const Point& rhs) { return !operator> (lhs,rhs); }
inline bool operator>= (const Point& lhs, const Point& rhs) { return !operator< (lhs,rhs); }


#endif /* POINT_H_ */
