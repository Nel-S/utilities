#ifndef _MATH_H
#define _MATH_H

extern const double U_PI, U_SQRT_2, U_SQRT_3, U_TWO_PI;

#ifndef max(X, Y)
    #define max(X,Y) ((X) > (Y) ? (X) : (Y))
#endif
#ifndef min(X, Y)
    #define min(X,Y) ((X) < (Y) ? (X) : (Y))
#endif
#ifndef IF_EXISTS(EXPRESSION_TO_TEST, OUTPUT_IF_EXPRESSION_ISNT_NULL)
    #define IF_EXISTS(EXPRESSION_TO_TEST, OUTPUT_IF_EXPRESSION_ISNT_NULL) ((EXPRESSION_TO_TEST) ? (OUTPUT_IF_EXPRESSION_ISNT_NULL) : 0)
#endif

#endif