/*
 * lakemaggiore.h
 *
 *  Created on: 20/may/2016
 *      Author: MatteoG
 */

#ifndef lakemaggiore_H_
#define lakemaggiore_H_

#include <math.h>
#include <vector>
#include <string>
#include "lake.h"

namespace std{

class lakemaggiore : public lake {
public:
        lakemaggiore();
        virtual ~lakemaggiore();

        /**
          * Level-Surface-Storage functions
          */

        double levelToSurface(double h);
        double levelToStorage(double h);
        double storageToLevel(double s);


protected:


        /**
          * Function to compute min-max release
          *  - s = lake storage
          *  - cday = day of the year (for MEF)
          */

        double min_release( double s, int cday );
        double max_release( double s, int cday );
        double s_max(int cday);

        vector<double> theta,theta1,theta2;
        vector<double> smax;



};

}

#endif
