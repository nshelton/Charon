out vec4 o_color;

uniform vec4 u_rotation;
uniform vec4 u_position;

uniform vec4 u_paramA;
uniform vec4 u_paramB;
uniform vec4 u_paramC;
uniform vec4 u_paramD;
uniform vec4 u_metaParam;
uniform vec4 u_metaParam2;

#define _LEVELS int(u_metaParam2.x * 10)

float thresh(float d) {
    return d * u_metaParam.x / 100;
}

float saturate(float x) {
    return clamp(x, 0.0, 1.0);
}

vec3 viridis_quintic(float x) {
    x = saturate(x);
    vec4 x1 = vec4(1.0, x, x * x, x * x * x); // 1 x x2 x3
    vec4 x2 = x1 * x1.w * x; // x4 x5 x6 x7
    return vec3(dot(x1.xyzw, vec4(+0.280268003, -0.143510503, +2.225793877, -14.815088879)) + dot(x2.xy, vec2(+25.212752309, -11.772589584)), dot(x1.xyzw, vec4(-0.002117546, +1.617109353, -1.909305070, +2.701152864)) + dot(x2.xy, vec2(-1.685288385, +0.178738871)), dot(x1.xyzw, vec4(+0.300805501, +2.614650302, -12.019139090, +28.933559110)) + dot(x2.xy, vec2(-33.491294770, +13.762053843)));
}

vec3 inferno_quintic(float x) {
    x = saturate(x);
    vec4 x1 = vec4(1.0, x, x * x, x * x * x); // 1 x x2 x3
    vec4 x2 = x1 * x1.w * x; // x4 x5 x6 x7
    return vec3(dot(x1.xyzw, vec4(-0.027780558, +1.228188385, +0.278906882, +3.892783760)) + dot(x2.xy, vec2(-8.490712758, +4.069046086)), dot(x1.xyzw, vec4(+0.014065206, +0.015360518, +1.605395918, -4.821108251)) + dot(x2.xy, vec2(+8.389314011, -4.193858954)), dot(x1.xyzw, vec4(-0.019628385, +3.122510347, -5.893222355, +2.798380308)) + dot(x2.xy, vec2(-3.608884658, +4.324996022)));
}

vec3 magma_quintic(float x) {
    x = saturate(x);
    vec4 x1 = vec4(1.0, x, x * x, x * x * x); // 1 x x2 x3
    vec4 x2 = x1 * x1.w * x; // x4 x5 x6 x7
    return vec3(dot(x1.xyzw, vec4(-0.023226960, +1.087154378, -0.109964741, +6.333665763)) + dot(x2.xy, vec2(-11.640596589, +5.337625354)), dot(x1.xyzw, vec4(+0.010680993, +0.176613780, +1.638227448, -6.743522237)) + dot(x2.xy, vec2(+11.426396979, -5.523236379)), dot(x1.xyzw, vec4(-0.008260782, +2.244286052, +3.005587601, -24.279769818)) + dot(x2.xy, vec2(+32.484310068, -12.688259703)));
}

vec3 plasma_quintic(float x) {
    x = saturate(x);
    vec4 x1 = vec4(1.0, x, x * x, x * x * x); // 1 x x2 x3
    vec4 x2 = x1 * x1.w * x; // x4 x5 x6 x7
    return vec3(dot(x1.xyzw, vec4(+0.063861086, +1.992659096, -1.023901152, -0.490832805)) + dot(x2.xy, vec2(+1.308442123, -0.914547012)), dot(x1.xyzw, vec4(+0.049718590, -0.791144343, +2.892305078, +0.811726816)) + dot(x2.xy, vec2(-4.686502417, +2.717794514)), dot(x1.xyzw, vec4(+0.513275779, +1.580255060, -5.164414457, +4.559573646)) + dot(x2.xy, vec2(-1.916810682, +0.570638854)));
}

vec2 hartverdrahtet(vec3 p) {

    vec3 cs = u_paramA.xyz;
    float fs = u_paramA.w;

    vec3 fc = u_paramB.xyz;
    float fu = u_paramB.w;
    float orbit = 0.0;

    float fd = 0.763;
    float dEfactor = 1.;
   //int fractal_iterations = 12;
    for(int i = 0; i < _LEVELS; i++) {
        vec3 start = p;
      //box folding
        p = 2. * clamp(p, -cs, cs) - p;
      //inversion
        float k = max(fs / dot(p, p), 1.);
        p *= k;
        dEfactor *= k;
      //julia seed
        p += fc;

        orbit += length(start - p);
    }
   //call basic shape and scale its DE
   //need to adjust fractal_distancemult with non zero julia seed
    float rxy = length(p.xy) - fu;

   //distance from pos to the pseudo kleinian basic shape ...

    return vec2(fd * max(rxy, abs(length(p.xy) * p.z) / sqrt(dot(p, p))) / abs(dEfactor), orbit);
}

// greetz 2 Mikael Hvidtfeldt Christensen
// http://blog.hvidtfeldts.net/index.php/2011/09/distance-estimated-3d-fractals-v-the-mandelbulb-different-de-approximations/
vec2 mandelbulb(vec3 pos) {

    float Power = u_paramA.x;
    float Bailout = u_paramA.y;

    vec3 z = pos; //
    float dr = 1.0;
    float r = 0.0;
    float iter;

    for(int i = 0; i < _LEVELS; i++) {
        r = length(z);
        if(r > Bailout)
            break;

        // convert to polar coordinates
        float theta = acos(z.z / r) + u_paramA.z;
        float phi = atan(z.y, z.x) + u_paramA.w;
        dr = pow(r, Power - 1.0) * Power * dr + 1.0;

        // scale and rotate the point
        float zr = pow(r, Power);
        theta = theta * Power;
        phi = phi * Power;

        // convert back to cartesian coordinates
        z = zr * vec3(sin(theta) * cos(phi), sin(phi) * sin(theta), cos(theta));
        z += pos;
        iter++;
    }

    return vec2(0.5 * log(r) * r / dr, iter);
}

vec2 tglad_variant(vec3 z0) {
    z0 = z0 / u_metaParam2.y;
    // z0 = modc(z0, 2.0);
    float mr = u_paramD.x, mxr = u_paramD.y;

    vec4 scale = vec4(u_paramD.z);
    vec4 p0 = u_paramA.xyzz;
    vec4 z = vec4(z0, 1.0);

    float orbit = 0;
    for(int n = 0; n < _LEVELS; n++) {
        vec3 start = z.xyz;
        z.xyz = clamp(z.xyz, -u_paramB.x, u_paramB.x) * 2.0 - z.xyz;
        z *= scale / clamp(dot(z.xyz, z.xyz), mr, mxr);
        z += p0;
        orbit += length(start - z.xyz);
    }
    float dS = (length(max(abs(z.xyz) - u_paramC.xyz * 2.0, 0.0)) - 0) / z.w;
    return vec2(u_metaParam2.y * dS, orbit);
}

//https://www.shadertoy.com/view/4ds3zn
vec2 appolonian(vec3 p) {
    float s = u_paramA.y * 2;
    float scale = u_paramA.x * 2;

    vec4 orb = vec4(1000.0);

    for(int i = 0; i < _LEVELS; i++) {
        p = -1.0 + 2.0 * fract(0.5 * p + 0.5);
        //p = rotate_vector(p, fromtwovectors(vec3(0,1,0), u_paramB.xyz));

        float r2 = dot(p, p);

        orb = min(orb, vec4(abs(p), r2));

        float k = s / r2;
        p *= k;
        scale *= k;
    }

    return vec2(0.25 * length(p) / scale, length(orb.y));
}

vec2 DE(in vec3 p) {

    float type = u_metaParam.a * 5.0;
    if(type < 1.0) {
        return appolonian(p);

    } else if(type < 2.0) {
        return tglad_variant(p);

    } else if(type < 3.0) {
        return hartverdrahtet(p);

    } else if(type < 4.0) {
        return mandelbulb(p);

    }
    return vec2(0);

}

vec3 calculate_normal(in vec3 p, in vec3 ro) {
    float e = thresh(length(p - ro)) / 2.0;
    const vec3 small_step = vec3(e, 0.0, 0.0);

    float gradient_x = DE(p + small_step.xyy).x - DE(p - small_step.xyy).x;
    float gradient_y = DE(p + small_step.yxy).x - DE(p - small_step.yxy).x;
    float gradient_z = DE(p + small_step.yyx).x - DE(p - small_step.yyx).x;

    vec3 normal = vec3(gradient_x, gradient_y, gradient_z);

    return normalize(normal);
}

vec3 getPalette(float t) {
    if(u_metaParam.y < 0.25) {
        return viridis_quintic(t);
    } else if(u_metaParam.y < 0.5) {
        return inferno_quintic(t);
    } else if(u_metaParam.y < 0.75) {
        return magma_quintic(t);
    } else {
        return plasma_quintic(t);
    }
}

vec3 ray_march(in vec3 ro, in vec3 rd) {
    float total_distance_traveled = 0.0;
    const int NUMBER_OF_STEPS = 100;
    const float MINIMUM_HIT_DISTANCE = 0.001;
    const float MAXIMUM_TRACE_DISTANCE = 1000.0;

    for(int i = 0; i < NUMBER_OF_STEPS; ++i) {
        vec3 current_position = ro + total_distance_traveled * rd;

        vec2 dist = DE(current_position);

        if(dist.x < thresh(total_distance_traveled)) {
            vec3 normal = calculate_normal(current_position, ro);
            vec3 light_position = vec3(2.0, -5.0, 3.0);
            vec3 direction_to_light = normalize(ro - current_position);

            float lambert = 0.5 * dot(normal, direction_to_light) + 0.5;
            // return vec3(1.0 - float(i) / float(NUMBER_OF_STEPS)) * 
                // inferno_quintic(dist.y * u_metaParam.y * 100.0) * diffuse_intensity;
            float t = 1.0 - float(i) / float(NUMBER_OF_STEPS);
            vec3 col = getPalette(t);
            return col * t;//* pow(lambert, 0.5);
            // return vec3( dist.y * u_metaParam.y);
        }

        if(total_distance_traveled > MAXIMUM_TRACE_DISTANCE) {
            break;
        }
        total_distance_traveled += dist.x * u_metaParam.z;
    }
    return vec3(0.0);
}

float rand(vec2 co) {
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}
// Quaternion multiplication
// http://mathworld.wolfram.com/Quaternion.html
vec4 qmul(vec4 q1, vec4 q2) {
    return vec4(q2.xyz * q1.w + q1.xyz * q2.w + cross(q1.xyz, q2.xyz), q1.w * q2.w - dot(q1.xyz, q2.xyz));
}

// Vector rotation with a quaternion
// http://mathworld.wolfram.com/Quaternion.html
vec3 rotate_vector(vec3 v, vec4 r) {
    vec4 r_c = r * vec4(-1, -1, -1, 1);
    return qmul(r, qmul(vec4(v, 0), r_c)).xyz;
}

void main() {
    vec2 uv = vUV.st * 2.0 - 1.0;
    uv.y *= uTDOutputInfo.res.x / uTDOutputInfo.res.y;
    uv *= 2.0;

    uv += vec2(0.5 + 0.5 * rand(uv + u_metaParam2.z) * 0.005, 0.5 + 0.5 * rand(uv + u_metaParam2.z + 131.7) * 0.005);

    vec3 ro = u_position.xyz;
    vec3 rd = vec3(uv, 1.0);

    rd = rotate_vector(rd, u_rotation);

    vec3 shaded_color = ray_march(ro, rd);

    o_color = vec4(shaded_color, 1.0);
}
