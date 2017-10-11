//
// Created by sachetto on 01/10/17.
//

#include "utils.h"

int binary_search(double **a, double num, int first, int last, int column) {
    int middle = first + ((last-first)/2);
    if (a[middle][column] == num)
        return middle;
    else if (first == last || first > last)
        return -1;
    else if (a[middle][column] < num)
        return binary_search(a,num,middle+1,last,column);
    else if (a[middle][column] > num)
        return binary_search(a,num,first,middle-1,column);

    return -1;
}

int inside_mesh(double **a, double x, double y, double z, int first, int last) {
    int xval = binary_search(a,x,first,last,0);
    if (xval == -1)
        return -1;
    else {
        double xnum = a[xval][0];
        int i = xval;
        while (a[i][0] == xnum && i != first)
            i--;
        int xfirst;
        if (i == first)
            xfirst = i;
        else
            xfirst = i + 1;
        int j = xval;
        while (a[j][0] == xnum && j != last)
            j++;
        int xlast;
        if (j == last)
            xlast = j;
        else
            xlast = j - 1;
        int yval = binary_search(a,y,xfirst,xlast,1);
        if (yval == -1)
            return -1;
        else {
            double ynum = a[yval][1];
            int k = yval;
            while (a[k][1] == ynum && k != first)
                k--;
            int yfirst;
            if (k == first)
                yfirst = k;
            else
                yfirst = k + 1;
            int l = yval;
            while (a[l][1] == ynum && l != last)
                l++;
            int ylast;
            if (l == last)
                ylast = l;
            else
                ylast = l - 1;
            int zval = binary_search(a,z,yfirst,ylast,2);
            if (zval == -1)
                return -1;
            else
                return zval;
        }

    }
}