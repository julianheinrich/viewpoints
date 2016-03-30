
varying vec3 packed_data_0 ;
varying vec4 packed_data_1 ;
varying vec4 packed_data_2 ;
varying vec4 packed_data_3 ;
varying vec4 packed_data_4 ;

//varying vec3 N;
#define COLOR packed_data_3
#define NORMAL normalize(packed_data_0.xyz)

uniform float fog_enabled;
varying float fog;

uniform sampler2D bgTextureMap;
uniform vec3 fogSolidColor;
uniform float fogIsSolidColor;

varying vec2 bgTextureLookup;

uniform bool lighting_enabled;
uniform bool two_sided_lighting_enabled;
uniform int light_count;
uniform vec4 interior_color;
uniform float interior_color_threshold;
uniform float shininess;
uniform float shininess_0;
uniform bool use_interior_color_threshold;
uniform int spec_count;
uniform float spec_value;
uniform float spec_value_0;

#include ANAGLYPH_HEADER

uniform float isStretched;
uniform float isCentered;
uniform float isCenteredOrRepeated;
uniform float isTiled;
uniform vec2 tileSize;
uniform vec2 tiledSize;
uniform vec2 viewImageSize;
uniform vec2 pixelSize;
uniform vec2 halfPixel;

// http://aras-p.info/blog/2009/07/30/encoding-floats-to-rgba-the-final/
vec4 EncodeFloatRGBA( float v ) {
  vec4 enc = vec4(1.0, 255.0, 65025.0, 160581375.0) * v;
  enc = fract(enc);
  enc -= enc.yzww * vec4(1.0/255.0,1.0/255.0,1.0/255.0,0.0);
  return enc;
}

void main()
{
  gl_FragColor = EncodeFloatRGBA(gl_FragCoord.z);
}

