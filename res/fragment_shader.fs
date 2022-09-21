#version 330

in vec2 fragTexCoord;
in vec4 fragColor;
in vec3 fragNormal;
in vec3 fragPosition;

flat in int drawGrid;

out vec4 finalColor;

const float offset = 1.0 / 20.0;
const float one_minus_offset = 1.0 - offset;
const vec4 pink = vec4(1, 0.46, 1, 1.0);
const float ambientStrength = 0.1;

struct Light {
    int enabled;
    int type;
    vec3 position;
    vec3 target;
    vec4 color;
};

// Input lighting values
uniform Light lights[1];

void main() {

    vec3 lightDot = vec3(0.0);
    vec3 normal = normalize(fragNormal);
    vec3 ambient = vec3(1.0);

    vec3 lp = lights[0].position;
    vec3 lt = lights[0].target;
    vec4 lc = lights[0].color;

    vec3 light = vec3(0.0);
    light = -normalize(lt - lp);
    float NdotL = max(dot(normal, light), 0.0);
    lightDot += lc.rgb * NdotL;
    ambient = ambientStrength * lc.rgb;

    bool draw_countour = fragTexCoord.x < offset || fragTexCoord.y < offset || fragTexCoord.x > one_minus_offset || fragTexCoord.y > one_minus_offset;
    bool selected = (fragColor.a == 0);

    finalColor = fragColor;

    if(drawGrid == 0) {
        if(selected) {
            if(draw_countour) {
                finalColor = pink;
            } else {
                finalColor = vec4(fragColor.rgb, 1);
            }
        } else {
            finalColor = fragColor;
        }
    } else if(drawGrid == 1) {
        if(draw_countour) {
            if(selected) {
                finalColor = pink;
            } else {
                finalColor = vec4(0,0,0,1);
            }
        }
        else {
            if(selected) {
                finalColor = vec4(fragColor.rgb, 1);
            } else {
                finalColor = fragColor;
            }
        }
    } else if(drawGrid == 2) {
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
    }

    vec3 result = (ambient + lightDot) * finalColor.rgb;
    finalColor = vec4(result, finalColor.a);
}


