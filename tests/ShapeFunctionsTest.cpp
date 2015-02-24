#include "catch.hpp"
extern "C" {
    #define restrict
    #include "interpolate.h"
    #undef restrict
}

TEST_CASE("Linear Shape Functions", "[linear][short]") {
    double phi[4];
    SECTION("Outside of Element") {
        SECTION("-X") {
            const double x_local = -1;
            const double y_local = 0;
            tent(phi, phi+1, phi+2, phi+3, x_local, y_local);
        }
        SECTION("+X") {
            const double x_local = 2;
            const double y_local = 0;
            tent(phi, phi+1, phi+2, phi+3, x_local, y_local);
        }
        SECTION("-Y") {
            const double x_local = 0;
            const double y_local = -1;
            tent(phi, phi+1, phi+2, phi+3, x_local, y_local);
        }
        SECTION("+Y") {
            const double x_local = 0;
            const double y_local = 2;
            tent(phi, phi+1, phi+2, phi+3, x_local, y_local);
        }
        CHECK(phi[0] == 0);
        CHECK(phi[1] == 0);
        CHECK(phi[2] == 0);
        CHECK(phi[3] == 0);
    }
}
