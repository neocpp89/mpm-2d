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
    EXIT_ERROR_CS_SOL,
    EXIT_ERROR_CS_DUP,

    EXIT_ERROR_CFG_PARSE=0x20,
    EXIT_ERROR_MATERIAL_FILE,
    EXIT_ERROR_BC_PROPS,

    EXIT_ERROR_DT_TOO_SMALL=0x40,

    EXIT_ERROR_THREADING=0x80
};

#endif /* __EXITCODES_H__ */
