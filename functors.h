#ifndef FUNCTORS_H
#define FUNCTORS_H

struct Compare_Points3i {
    bool operator() (const cv::Point3i& a, const cv::Point3i& b) const {
        if ( a.x < b.x )
            return true;
        if ( a.x == b.x && a.y < b.y)
            return true;
        if (a.x == b.x && a.y == b.y)
            return (a.z < b.z);

        return false;
    }
};

#endif // FUNCTORS_H
