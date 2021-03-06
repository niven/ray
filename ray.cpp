/*
clang++ -std=c++11 -g -Wall -Wextra -pedantic -Wpadded -Wno-missing-braces -Wno-gnu-anonymous-struct -O3 ray.cpp -o ray; and ./ray
-Wpadded? yes, to find errors in bmp struct packing
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define ARRAY_COUNT(array) ( sizeof(array) / sizeof(array[0]) )

#define f32 float
#define F32MAX FLT_MAX

#define u8 unsigned char
#define u16 unsigned short
#define u32 unsigned int
#define s8 signed char
#define s16 signed short
#define s32 signed int

#define internal static

typedef union v3 {
	struct {
		f32 x,y,z;
	};
	
	struct {
		f32 r,g,b;
	};

	struct {
		f32 u,v,w;
	};
	
	f32 E[3];
} v3;

internal inline v3 operator*( f32 f, v3 v ) {
	v3 result = { f* v.x, f*v.y, f*v.z };
	return result;
}
internal inline v3 operator/( v3 a, v3 b ) {
	v3 result = { a.x/b.x, a.y/b.y, a.z/b.z };
	return result;
}
internal inline v3 operator/( v3 a, f32 b ) {
	v3 result = { a.x/b, a.y/b, a.z/b };
	return result;
}
internal inline v3 operator-( v3 a, v3 b ) {
	v3 result = { a.x-b.x, a.y-b.y, a.z-b.z };
	return result;
}
internal inline v3 operator+( v3 a, v3 b ) {
	v3 result = { a.x+b.x, a.y+b.y, a.z+b.z };
	return result;
}
internal inline v3 operator+=( v3 &a, v3 b ) {
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	return a;
}

internal inline v3 operator-( v3 a ) {
	v3 result;
	result.x = -a.x;
	result.y = -a.y;
	result.z = -a.z;
	return result;
}

internal v3 V3( f32 x, f32 y, f32 z ) {
	v3 result = { x, y, z };
	return result;
}

inline internal f32 inner( v3 a, v3 b ) {
	f32 result = a.E[0] * b.E[0] + a.E[1] * b.E[1] + a.E[2] * b.E[2];
	return result;
}

internal v3 cross( v3 a, v3 b ) {
	v3 result;
	
	result.x = a.y*b.z - a.z*b.y;
	result.y = a.z*b.x - a.x*b.z;
	result.z = a.x*b.y - a.y*b.x;
	
	return result;
}

internal inline v3 NOZ( v3 n ) {
	v3 result;
	f32 d = sqrt( inner(n,n) );
	// What if d is 0??
	result.x = n.x / d;
	result.y = n.y / d;
	result.z = n.z / d;
	return result;
}

internal inline v3 hadamard( v3 a, v3 b ) {
	v3 result = { a.x*b.x, a.y*b.y, a.z*b.z };
	return result;
}

internal v3 lerp( v3 a, f32 d, v3 b ) {
	f32 c = 1.0f - d;
	v3 result = { c*a.x + d*b.x, c*a.y + d*b.y, c*a.z + d*b.z };
	return result;
}

// don't pad so we write a correct header!
#pragma pack(push, 1)
typedef struct bitmap_file_header {
	u16 header_field;
	u32 file_size;
	u16 reserved_0;
	u16 reserved_1;
	u32 image_offset;
	/* BITMAPINFOHEADER */
	u32 header_size;
	s32 image_width;
	s32 image_height;
	u16 color_planes;
	u16 bits_per_pixel; /* Typical values are 1, 4, 8, 16, 24 and 32. */
	u32 compression_method;
	u32 image_size;
	s32 hresolution; /* pix per m */
	s32 vresolution;
	u32 palette_colors;
	u32 important_colors;
} bitmap_file_header;
#pragma pack(pop)

typedef struct plane {
	v3 N; // normal
	f32 d; // offset from origin
	u32 material_index;
} plane;

typedef struct sphere {
	v3 P;
	f32 r;
	u32 material_index;
} sphere;

typedef struct ellipsoid {
	v3 P;
	v3 a;
	u32 material_index;
} ellipsoid;

typedef struct material {
	f32 scatter;
	v3 color_reflect;
	v3 color_emit;
} material;

typedef struct tri {
	union {
		v3 vert[3];
		struct {
			v3 a;
			v3 b;
			v3 c;
		};
	};

	v3 N; // which side is up??
	u32 material_index[3]; // 1 for each vertex
} tri;



typedef struct world {

	material *materials;

	plane *planes;
	sphere *spheres;
	ellipsoid *ellipsoids;
	tri *triangles;
	
	u32 material_count;
	u32 plane_count;
	u32 sphere_count;
	u32 ellipsoid_count;
	u32 triangle_count;
	
} world;

typedef struct image_rgba {
	u32 width;
	u32 height;
	u32 *pixels;
} image_rgba;

void write_bmp( const char* filename, image_rgba img ) {
	
	struct bitmap_file_header bmf;
	bmf.header_field = 0x4d << 8 | 0x42;
	bmf.reserved_0 = 0;
	bmf.reserved_1 = 0;
	bmf.image_offset = sizeof(bitmap_file_header);
	
	bmf.color_planes = 1;
	bmf.bits_per_pixel = 32; /* Typical values are 1, 4, 8, 16, 24 and 32. */
	bmf.header_size = sizeof(bitmap_file_header) - 14;
	bmf.image_width = img.width;
	bmf.image_height = img.height; // flip vertically so y increases down
	bmf.compression_method = 0;
	bmf.hresolution = 1000; /* pix per m */
	bmf.vresolution = 1000;
	bmf.palette_colors = 0;
	bmf.important_colors = 0;
	bmf.image_size = bmf.image_width * bmf.image_height * 4; // RGBA
	bmf.file_size = sizeof(bitmap_file_header) + bmf.image_size;
	
	FILE *bmp = fopen(filename, "wb");
	
	fwrite(&bmf, sizeof(bitmap_file_header), 1, bmp);
	fwrite(img.pixels, bmf.image_size, 1, bmp);
	
	fclose( bmp );
}

image_rgba allocate_image( u32 width, u32 height ) {
	image_rgba result = { width, height, 0 };
	result.pixels = (u32*)malloc( sizeof(u32) * width * height );
	return result;
}

// gamma correction, see wikipedia
internal f32 sRGB_from_linear( f32 linear ) {
	f32 S;
	
	if( linear < 0 ) {
		linear = 0;
	}
	if( linear > 1.0f ) {
		linear = 1.0f;
	}
	
	if( linear <= 0.0031308f ) {
		S = linear * 12.92f;
	} else {
		S = 1.055f * pow(linear, 1.0f/2.4f) - 0.055f;
	}
	return S;
}

f32 random_unilateral( void ) {
	f32 result = (f32)rand() / (f32) RAND_MAX;
	return result;
}

f32 random_bilateral( void ) {
	f32 result = -1.0f + 2.0f * random_unilateral();
	return result;
}


f32 ray_intersects_plane( v3 ray_origin, v3 ray_direction, plane p ) {
	f32 result;
	
	f32 tolerance = 0.000001f;
	f32 denom = inner( p.N, ray_direction );
	if( denom < -tolerance || denom > tolerance ) {
		result = ( -p.d  - inner( p.N, ray_origin ) ) / denom;		
	} else {
		result = F32MAX;
	}

	return result;
}


u8 point_in_triangle( v3 p, tri t, v3 &uvw ) {

	// Compute vectors        
	v3 v0 = t.c - t.a;
	v3 v1 = t.b - t.a;
	v3 v2 = p - t.a;

	// Compute dot products
	f32 dot00 = inner(v0, v0);
	f32 dot01 = inner(v0, v1);
	f32 dot02 = inner(v0, v2);
	f32 dot11 = inner(v1, v1);
	f32 dot12 = inner(v1, v2);

	// Compute barycentric coordinates
	f32 inv_denom = 1.0f / (dot00 * dot11 - dot01 * dot01);
	uvw.u = (dot11 * dot02 - dot01 * dot12) * inv_denom;
	uvw.v = (dot00 * dot12 - dot01 * dot02) * inv_denom;
	uvw.w = 1.0f - uvw.u - uvw.v;

	// Check if point is in triangle
	u8 result = (uvw.u >= 0) && (uvw.v >= 0) && (uvw.w >= 0);
	return result;
}

f32 ray_intersects_triangle( v3 ray_origin, v3 ray_direction, tri t, v3 &uvw ) {

	// ray intersect plane defined by t
	// subproblem: create normal from points of triangle which I think is:
	// cross(b-a, c-a); // make vert a the origin, then cross product b and c should be a normal
	
	// ray-plane:
	f32 d; // r_0 + d * r_d
	f32 denom = inner( t.N, ray_direction );
	f32 tolerance = 0.000001f;
	
	if( denom < -tolerance || denom > tolerance ) {
		// triangle center???
		// how do we get the offset from the origin? Do we even need it?
		d = ( - inner( t.N, ray_origin ) ) / denom;	
		
		v3 intersect_p = ray_origin + d * ray_direction;
		// check is within bounds of actual triangle
		u8 hit = point_in_triangle( intersect_p, t, uvw );
		if( hit ) {
			return d;
		} else {
			return F32MAX;
		}
		
	}
	
	return F32MAX;
}

f32 ray_intersects_sphere( v3 ray_origin, v3 ray_direction, sphere s ) {
	f32 result;

	/*
		point on a sphere: x^2 + y^2 + z^2 = r^2  ->  inner(p,p) - r^2 = 0
		point on a ray: r_0 + t*r_d
	
		translating the sphere to the origin is just moving -s.origin
	
		inner( t*r_d + (r_0 - s_0), t*r_d + (r_0 - s_0) ) - s.r*s.r = 0
		inner product is distributive
		inner(r_d,_rd) * t^2 + 2 * inner( (r_0 - s_0), r_d ) * t + inner( (r_0 - s_0), (r_0 - s_0) ) - s.r*s.r = 0
		(r_0 - s_0) is basically relative_origin
		Solving quadratic formula:
		x_12 = ( -b +/- sqrt( b^2 - 4*a*c ) ) / 2*a
		a = inner(r_d,_rd)
		b = 2 * inner( relative_origin, r_d )
		c = inner( relative_origin, relative_origin ) - s.r*s.r
	*/
	
	v3 relative_origin = ray_origin - s.P;
	// TODO: maybe think about this more. ray_direction is always normalized, but it need not be
	// mabe normalizing it always is more expensive than saving this is worth
	f32 a = 1.0f; // inner( ray_direction, ray_direction );
	f32 b = 2.0f * inner( relative_origin, ray_direction );
	f32 c = inner( relative_origin, relative_origin ) - s.r*s.r;
	
	// we divide by 2a which is only 0 is ray_direction is 0, so that won't happen
	// if we would sqrt() a negative it means no solution meaning no intersection
	f32 under = b*b - 4*a*c;
	f32 tolerance = 0.000001f;
	f32 a2 = 2.0f * a;
	if( under < 0 || (a2 < tolerance && a2 > -tolerance) ) {
		return F32MAX;
	}
	
	f32 root = sqrt(under);
	f32 x1 = ( -b - root ) / a2;
	f32 x2 = ( -b + root ) / a2;
	
	f32 min_hit_distance = 0.0001f;
	
	if( x1 > min_hit_distance && x1 < x2 ) {
		result = x1;
	} else if( x2 > min_hit_distance) {
		result = x2;
	} else {
		result = F32MAX;
	}

	return result;
}

f32 ray_intersects_ellipsoid( v3 ray_origin, v3 ray_direction, ellipsoid e ) {
	f32 result;
	v3 relative_origin = ray_origin - e.P;
	
	/*
		point on an ellipsoid: x^2/a0^2 + y^2/a1^2 + z^2/a2^2 = 1  ->  inner( p/a, p/a ) - 1 = 0
		point on a ray: r_0 + t*r_d

		translating the ellipsoid to the origin is just moving -e.origin
		inner( (r_0 + t*r_d)/a , (r_0 + t*r_d)/a) - 1 = 0
		inner( r_d/a, r_d/a ) * t^2 + 2 * inner( r_d/a, (r_0 - s_0)/a ) * t + inner( (r_0 - s_0)/a, (r_0 - s_0)/a ) - 1 = 0
		a = inner( r_d/a, r_d/a )
		b = 2 * inner( r_d/a, (r_0 - s_0)/a )
		c = inner( (r_0 - s_0)/a, (r_0 - s_0)/a ) - 1
	*/

	v3 rd_over_a = ray_direction / e.a;
	v3 ro_over_a = relative_origin / e.a;

	f32 a = inner( rd_over_a, rd_over_a );
	f32 b = 2 * inner( rd_over_a, ro_over_a );
	f32 c = inner( ro_over_a, ro_over_a ) - 1;
	
	// we divide by 2a which is only 0 is ray_direction is 0, so that won't happen
	// if we would sqrt() a negative it means no solution meaning no intersection
	f32 under = b*b - 4*a*c;
	f32 tolerance = 0.000001f;
	f32 a2 = 2.0f * a;
	if( under < 0 || (a2 < tolerance && a2 > -tolerance) ) {
		return F32MAX;
	}
	
	f32 root = sqrt(under);
	f32 x1 = ( -b - root ) / a2;
	f32 x2 = ( -b + root ) / a2;
	
	f32 min_hit_distance = 0.0001f;
	
	if( x1 > min_hit_distance && x1 < x2 ) {
		result = x1;
	} else if( x2 > min_hit_distance) {
		result = x2;
	} else {
		result = F32MAX;
	}

	return result;
}

v3 ray_cast( world* w, v3 ray_origin, v3 ray_direction ) {

	v3 result = {};
	
	f32 hit_distance;
	f32 min_hit_distance = 0.0001f;
	u8 hit;
	material mat_value = w->materials[0];
	
	u32 bounces_per_ray = 64;
	
	f32 t;
	v3 attenuation = V3(1,1,1);
	for( u32 bounce_count=0; bounce_count<bounces_per_ray; bounce_count++ ) {
		
		hit_distance = F32MAX;
		hit = 0;
		mat_value = w->materials[0];
		v3 next_normal;
		
		for( u32 plane_index=0; plane_index < w->plane_count; plane_index++ ) {
			plane current_plane = w->planes[plane_index];
		
			t = ray_intersects_plane( ray_origin, ray_direction, current_plane );
			if( t > min_hit_distance && t < hit_distance ) {
				hit_distance = t;
				mat_value = w->materials[current_plane.material_index];
				hit = 1;
				next_normal = current_plane.N;
			}
		}

		for( u32 sphere_index=0; sphere_index < w->sphere_count; sphere_index++ ) {
			sphere current_sphere = w->spheres[sphere_index];
		
			t = ray_intersects_sphere( ray_origin, ray_direction, current_sphere );
			if( t > min_hit_distance && t < hit_distance ) {
				hit_distance = t;
				mat_value = w->materials[current_sphere.material_index];
				hit = 1;
				next_normal = NOZ( ray_origin + t*ray_direction - current_sphere.P );
			}
		}

		for( u32 ellipsoid_index=0; ellipsoid_index < w->ellipsoid_count; ellipsoid_index++ ) {
			ellipsoid current_ellipsoid = w->ellipsoids[ellipsoid_index];
		
			t = ray_intersects_ellipsoid( ray_origin, ray_direction, current_ellipsoid );
			if( t > min_hit_distance && t < hit_distance ) {
				hit_distance = t;
				mat_value = w->materials[current_ellipsoid.material_index];
				hit = 1;
				// The normal of an ellipsoid is the gradient of x^2/a^ + y^2/b^2 + z^2/c^2 - 1 = 0
				// Which is (2x/a^2, 2y/b^2, 2z/c^2)
				// = 2 * point_on_surface / hadamard(e.a, e.a) (we can drop the 2 since we normalize anyway)
				// for a sphere that is just a vector from the center to the point on the surface and had( (1,1,1), (1,1,1) ) is (1,1,1)...
				// so there we translate the center of the sphere to the origin
				// for an ellipsoid we also need to translate something?
				// YES! We need to translate the ellipsoid to the origin, so the point on the surface will have the correct sign
				// when we multiply with the 1/a So actually moving the point is what we need to do, not the resulting normal
				v3 hit_location = ray_origin + t*ray_direction;
				next_normal = NOZ( (hit_location - current_ellipsoid.P) / hadamard(current_ellipsoid.a, current_ellipsoid.a) );
			}
		}
		
		for( u32 triangle_index=0; triangle_index < w->triangle_count; triangle_index++ ) {
			tri current_triangle = w->triangles[triangle_index];

			v3 uvw = {};
			t = ray_intersects_triangle( ray_origin, ray_direction, current_triangle, uvw );
			if( t > min_hit_distance && t < hit_distance ) {
				hit_distance = t;
				hit = 1;

				// f32 scatter_avg =
				// 	uwv.u * w->materials[current_triangle.material_index[0]].scatter +
				// 	uwv.v * w->materials[current_triangle.material_index[1]].scatter +
				// 	uwv.w * w->materials[current_triangle.material_index[2]].scatter;
				// scatter_avg /= 3.0f;
				//
				material m0 = w->materials[current_triangle.material_index[0]];
				material m1 = w->materials[current_triangle.material_index[1]];
				material m2 = w->materials[current_triangle.material_index[2]];
				v3 color_emit_avg =
					uvw.u * m0.color_emit +
					uvw.v * m1.color_emit +
					uvw.w * m2.color_emit;
				// mat_value.scatter = scatter_avg;
				mat_value.color_emit = NOZ(color_emit_avg);
				v3 color_reflect_avg =
					uvw.u * m0.color_reflect +
					uvw.v * m1.color_reflect +
					uvw.w * m2.color_reflect;
				mat_value.color_reflect = NOZ(color_reflect_avg);
				// printf("%f %f %f -> ce: %f %f %f\n", uvw.u, uvw.v, uvw.w, color_emit_avg.r, color_emit_avg.g, color_emit_avg.b);
				// mat_value = w->materials[3];
				next_normal = NOZ( current_triangle.N ); // NOZ( ray_origin + t*ray_direction - current_sphere.P );
			}
		}
		
		result += hadamard( attenuation, mat_value.color_emit );
		if( hit ) {
			// f32 cos_attenuation = inner( -ray_direction, next_normal );
			// if( cos_attenuation < 0 ) {
			// 	cos_attenuation = 0;
			// }
			attenuation = hadamard( attenuation, mat_value.color_reflect );

			ray_origin += hit_distance * ray_direction;

			v3 bounce_pure = ray_direction - 2.0f * inner(ray_direction, next_normal) *  next_normal;
			v3 bounce_random = NOZ( next_normal + V3(random_bilateral(), random_bilateral(), random_bilateral() ));
			ray_direction = NOZ( lerp( bounce_random, mat_value.scatter, bounce_pure ) );
		} else {
			break; // no more hits
		}
	}
	
	
	return result;
}

u32 rgbapack_4x8( v3 color ) {
	u32 result = 0xff << 24;
	result |= (u8)(0xff * color.r) << 16;
	result |= (u8)(0xff * color.g) << 8;
	result |= (u8)(0xff * color.b);
	return result;
}

u32* get_pixel_pointer(image_rgba img, u32 x, u32 y) {
	u32 *result = img.pixels + x + y * img.width;
	return result;
}

void raytrace_tile( world *w, image_rgba img, u32 x_start, u32 x_max, u32 y_start, u32 y_max ) {

	u32 rays_per_pixel = 32;
	v3 camera_p = V3(0, -10, 2.4f);
	// v3 camera_p = V3(0, 0.1f, 8);

	// look at the origin
	v3 camera_z = NOZ(camera_p);
	v3 camera_x = NOZ(cross( V3(0,0,1), camera_z ));
	v3 camera_y = NOZ(cross( camera_z, camera_x ));
	
	f32 film_dist = 1.0f;
	f32 film_w = 1.0f;
	f32 film_h = 1.0f;
	if( img.width > img.height ) {
		film_h = film_w * ( (f32)img.height / (f32)img.width);
	} else if( img.height > img.width ) {
		film_w = film_h * ( (f32)img.width / (f32)img.height);
	}
	
	f32 pixel_half_width = 0.5f / (f32)img.width;
	f32 pixel_half_height = 0.5f / (f32)img.height;
	f32 half_film_w = 0.5f * film_w;
	f32 half_film_h = 0.5f * film_h;
	// convert these from image space to film space which is a square from -1 to 1
	v3 film_center = camera_p - film_dist * camera_z;


	for( u32 y=y_start; y<y_max; y++ ) {
		f32 film_y = -1.0f + 2.0f * ( (f32)y / (f32)img.height );
		u32* out = get_pixel_pointer(img, x_start, y);
		for( u32 x=x_start; x<x_max; x++ ) {

			f32 film_x = -1.0f + 2.0f * ( (f32)x / (f32)img.width );

			v3 color = {};
			f32 contribution = 1.0f / (f32)rays_per_pixel;
			for( u32 ray_index = 0; ray_index<rays_per_pixel; ray_index++ ) {

				f32 offset_x = film_x + random_bilateral() * pixel_half_width;	
				f32 offset_y = film_y + random_bilateral() * pixel_half_height;	
				v3 film_p = film_center + offset_x * half_film_w * camera_x + offset_y * half_film_h * camera_y;
				v3 ray_origin = camera_p;
				v3 ray_direction = NOZ( film_p - camera_p );
				color += contribution * ray_cast( w, ray_origin, ray_direction );
			}
			
			v3 color_gamma_corrected = {
				sRGB_from_linear( color.r ),
				sRGB_from_linear( color.g ),
				sRGB_from_linear( color.b ),
			};
			u32 color_bmp = rgbapack_4x8( color_gamma_corrected );
			*out++ = color_bmp;
		}
	}
	
}


int main() {

	setbuf(stdout, NULL); // turn off printf buffering

	// image_rgba img = allocate_image( 2560, 1600 );
	image_rgba img = allocate_image( 1280, 720 );
// #if 0
	material materials[9] = {};
	materials[0].color_emit = V3(0.4f, 0.5f, 0.6f);
	// materials[0].color_reflect = V3(0.9f, 0.9f, 0.9f);
	materials[1].color_reflect = V3(.7f,.7f,.7f);
	materials[1].scatter = 0.2f;
	materials[2].color_reflect = V3(0.7f,0.5f,.3f);
	materials[2].color_emit = V3(0,0,0);

	// red light
	materials[3].color_reflect = V3(0.5f,0.5f,0.5f);
	materials[3].color_emit = V3(4.0f,0,0);
	// green light
	materials[4].color_reflect = V3(0.5f,0.5f,0.5f);
	materials[4].color_emit = V3(0,4.0f,0);
	// blue light
	materials[5].color_reflect = V3(0.5f,0.5f,0.5f);
	materials[5].color_emit = V3(0,0,4.0f);
	// orange light
	materials[6].color_reflect = V3(0.5f,0.5f,0.5f);
	materials[6].color_emit = V3(1.0f,0.4f,0.2f);
	// purple light
	materials[7].color_reflect = V3(0.5f,0.5f,0.5f);
	materials[7].color_emit = V3(100.0/255.0,0.0f,181.0/255.0);
	// mirror
	materials[8].color_reflect = V3(0.7f, 1.0f, 0.7f);
	materials[8].scatter = 1.0f;
	
	plane planes[1];
	planes[0].N = V3( 0, 0, 1.0f );
	planes[0].d = 0;
	planes[0].material_index = 1;
	
	sphere spheres[0];
#if 0
	spheres[0].r = 0.5;
	spheres[0].P = V3(-1,-2,0);
	spheres[0].material_index = 5;
	spheres[1].r = 0.5;
	spheres[1].P = V3(1,-2,0);
	spheres[1].material_index = 6;
	spheres[1].r = 0.4;
	spheres[1].P = V3(3,-2,1);
	spheres[1].material_index = 3;
	spheres[2].r = 1;
	spheres[2].P = V3(-1,-1,2);
	spheres[2].material_index = 4;
	spheres[3].r = 1;
	spheres[3].P = V3(1,-1,3);
	spheres[3].material_index = 5;
	// coordspheres
	u32 sidx = 4;
	for( s32 x=-2; x<3; x++ ) {
		for( s32 y=-2; y<3; y++ ) {
			spheres[sidx].r = 0.1; spheres[sidx].P = V3(x,y,0); spheres[sidx].material_index = 6;
			sidx++;
		}
	}

	spheres[29].r = 0.2;
	spheres[29].P = V3(1,-2,.5f);
	spheres[29].material_index = 7;
#endif
	
	tri triangles[1];
	triangles[0].a = V3(-1.0f, 0, 0);
	triangles[0].b = V3( 1.0f, 0, 0);
	triangles[0].c = V3( 0, 0, 2);
	triangles[0].N = -NOZ( cross( triangles[0].b - triangles[0].a, triangles[0].c - triangles[0].a) );
	triangles[0].material_index[0] = 3;
	triangles[0].material_index[1] = 4;
	triangles[0].material_index[2] = 5;
#if 0	
	triangles[1].a = V3( 1.0f, 0, 0);
	triangles[1].b = V3( 3.0f, 0, 0);
	triangles[1].c = V3( 2.0f, 0, 2);
	triangles[1].N = NOZ( cross( triangles[0].b - triangles[0].a, triangles[0].c - triangles[0].a) );
	triangles[1].material_index[0] = 3;
	triangles[1].material_index[1] = 4;
	triangles[1].material_index[2] = 5;
#endif
	
	
	// spheres[30].r = 1;
	// spheres[30].P = V3(1,-.5,1);
	// spheres[30].material_index = 5;
	
	
	// TODO: needs a direction / rotation!
	ellipsoid ellipsoids[0];
#if 0	
	ellipsoids[0].a = V3(2,1,1);
	ellipsoids[0].P = V3(1,-.5,1);
	ellipsoids[0].material_index = 8;
#endif
	
// #endif
#if 0
	material materials[6] = {};
	materials[0].color_emit = V3(1,1,1);

	materials[1].color_reflect = V3(0.6f, 0.6f, 0.6f);

	materials[2].color_emit = V3(5,0,0);
	materials[3].color_emit = V3(0,5,0);
	materials[4].color_emit = V3(0,0,5);
		
	plane planes[1];
	planes[0].N = V3( 0, 0, 1.0f );
	planes[0].d = 0;
	planes[0].material_index = 1;
	
	sphere spheres[3];
	spheres[0].r = 0.6;
	spheres[0].P = V3(-1,0,1);
	spheres[0].material_index = 2;
	spheres[1].r = 0.6;
	spheres[1].P = V3(.7f,.7f,1);
	spheres[1].material_index = 3;
	spheres[2].r = 0.6;
	spheres[2].P = V3(.7f,-.7f,1);
	spheres[2].material_index = 4;
#endif	
	
	
	world w;
	w.material_count = ARRAY_COUNT(materials);
	w.materials = materials;
	w.plane_count = ARRAY_COUNT(planes);
	w.planes = planes;
	w.sphere_count = ARRAY_COUNT(spheres);
	w.spheres = spheres;
	w.ellipsoid_count = ARRAY_COUNT(ellipsoids);
	w.ellipsoids = ellipsoids;
	w.triangle_count = ARRAY_COUNT(triangles);
	w.triangles = triangles;
	
	// let's do this!
	// so we pass a ray from every point on the image in the camera direction onto the world and see what we hit.
	// if we hit nothing, we use the mat[0] to set the pixel
	
	// raytrace_tile( &w, img, 0, img.width, 10, img.height );
	u8 core_count = 8;
	u32 tile_width = (img.width + core_count - 1) / core_count;
	u32 tile_height = tile_width;
	printf("Raytracing tiles %dx%d\n", tile_width, tile_height);

	for( u32 x_start = 0; x_start<img.width; x_start += tile_width ) {
		for( u32 y_start = 0; y_start<img.height; y_start += tile_height ) {
			u32 x_end = x_start + tile_width;
			u32 y_end = y_start + tile_height;
			if( x_end > img.width ) {
				x_end = img.width;
			}
			if( y_end > img.height ) {
				y_end = img.height;
			}
			printf("Raytracing tile [%d-%d,%d-%d]\n", x_start, x_end, y_start, y_end);
			raytrace_tile( &w, img, x_start, x_end, y_start, y_end );
		}
	}
	// raytrace_tile( &w, img, 0, img.width/2, 0, img.height/2 );
	// raytrace_tile( &w, img, img.width/2, img.width, 0, img.height/2 );
	// raytrace_tile( &w, img, 0, img.width/2, img.height/2, img.height );
	// raytrace_tile( &w, img, img.width/2, img.width, img.height/2, img.height );
	
	write_bmp( "test.bmp", img );
	
	printf("\nDone.\n");
	return 0;
}

