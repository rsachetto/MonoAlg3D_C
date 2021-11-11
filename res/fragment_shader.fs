#version 330

in vec2 fragTexCoord;
in vec4 fragColor;
flat in int drawGrid;
out vec4 finalColor;

const float offset = 1.0 / 20.0;
const float one_minus_offset = 1.0 - offset;

void main() {

    bool draw_countour = fragTexCoord.x < offset || fragTexCoord.y < offset || fragTexCoord.x > one_minus_offset || fragTexCoord.y > one_minus_offset;

    vec4 pink = vec4(1, 0.46, 1, 1.0);

    bool selected = (fragColor.a == 0);

    if(drawGrid == 0) {

        if(selected) {
            if(draw_countour) {
                finalColor = pink;
            }

            else {
                finalColor = vec4(fragColor.rgb, 1);
            }
        }

        else {
            finalColor = fragColor;
        }

        return;
    }

    if(drawGrid == 1) {
        if(draw_countour) {
            if(selected) {
                finalColor = pink;
            }
            else {
                finalColor = vec4(0,0,0,1);
            }
        }
        else {
            if(selected) {
                finalColor = vec4(fragColor.rgb, 1);
            }
            else {
                finalColor = fragColor;
            }
        }
        return;
    }

    if(drawGrid == 2) {
        if(draw_countour) {
            if(selected) {
                finalColor = pink;
            }
            else {
                finalColor = fragColor;
            }
        }
        else {
            finalColor = vec4(0,0,0,0);
        }
        return;
    }
}


