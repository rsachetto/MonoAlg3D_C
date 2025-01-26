#ifndef __GUI_COLORS_H
#define __GUI_COLORS_H

#include "../3dparty/raylib/src/raylib.h"
#define NUM_SCALES 7
#define NUM_COLORS 257

Color get_color(float value, int alpha, int current_scale);

#endif /* __GUI_COLORS_H */
