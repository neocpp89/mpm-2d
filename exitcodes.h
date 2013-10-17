/**
    \file exitcodes.h
    \author Sachith Dunatunga
    \date 16.10.13

    Contains an enum for the error exitcodes used by mpm_2d.
*/

#ifndef __EXITCODES_H__
#define __EXITCODES_H__

enum exitcodes_e {
    EXIT_ERROR_CS_ENTRY=0x01,
    EXIT_ERROR_CS_SOL=0x02,
};

#endif /* __EXITCODES_H__ */
