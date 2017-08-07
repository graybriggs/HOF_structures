#ifndef FUNCTORS_H
#define FUNCTORS_H

struct Compare_Points3i {
    bool operator() (cv::Point3i const& a, cv::Point3i const& b) const {
        if ( a.x < b.x ) return true;
        if ( a.x == b.x and a.y < b.y ) return true;
        if ( a.x == b.x and a.y == b.y ) return ( a.z < b.z );
        return false;
    }
};

#endif // FUNCTORS_H
