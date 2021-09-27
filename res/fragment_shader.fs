#version 330

in vec2 fragTexCoord;
in vec4 fragColor;
flat in int drawGrid;
out vec4 finalColor;

const float offset = 1.0 / 10.0;
const float one_minus_offset = 1.0 - offset;

void main() {

    if(drawGrid == 0) {
        finalColor = fragColor;
    }

    else {
        bool draw_countour = fragTexCoord.x < offset || fragTexCoord.y < offset || fragTexCoord.x > one_minus_offset || fragTexCoord.y > one_minus_offset;


        if(drawGrid == 1) {
            if(draw_countour) {
                finalColor = vec4(0,0,0,1);
            }
            else {
                finalColor = fragColor;
            }

        }

        if(drawGrid == 2) {
            if(draw_countour) {
                finalColor = fragColor;
            }
            else {
                finalColor = vec4(1,1,1,0);
            }

        }
    }
}


