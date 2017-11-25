/*
clang -std=c11 -g -Wpadded -Wall -Wextra -pedantic -Wno-missing-braces -Wno-gnu-anonymous-struct -O0 ray.c -o ray
-Wpadded? yes, to find errors in bmp struct packing
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

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
	
	f32 E[3];
} v3;

v3 operator*( f32 f, v3 v ) {
	v3 result = { f* v.x, f*v.y, f*v.z };
	return result;
}
v3 operator-( v3 a, v3 b ) {
	v3 result = { a.x-b.x, a.y-b.y, a.z-b.z };
	return result;
}
v3 operator+( v3 a, v3 b ) {
	v3 result = { a.x+b.x, a.y+b.y, a.z+b.z };
	return result;
}

internal v3 V3( f32 x, f32 y, f32 z ) {
	v3 result = { x, y, z };
	return result;
}

internal f32 inner( v3 a, v3 b ) {
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

internal v3 NOZ( v3 n ) {
	v3 result;
	f32 d = sqrt( inner(n,n) );
	// What if d is 0??
	result.x = n.x / d;
	result.y = n.y / d;
	result.z = n.z / d;
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

typedef struct material {
	v3 color;
} material;

typedef struct world {
	u32 material_count;
	material *materials;

	u32 plane_count;
	plane *planes;
	
	u32 sphere_count;
	sphere *spheres;
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
	f32 a = inner( ray_direction, ray_direction );
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
	
	f32 x1 = ( -b - sqrt(under) ) / a2;
	f32 x2 = ( -b + sqrt(under) ) / a2;
	
	if( x1 < x2 ) {
		result = x1;
	} else {
		result = x2;
	}

	return result;
}

v3 ray_cast( world* w, v3 ray_origin, v3 ray_direction ) {
	v3 result = w->materials[0].color;
	
	f32 hit_distance = F32MAX;
	
	for( u32 plane_index=0; plane_index < w->plane_count; plane_index++ ) {
		plane current_plane = w->planes[plane_index];
		
		f32 this_distance = ray_intersects_plane( ray_origin, ray_direction, current_plane );
		if( this_distance > 0 && this_distance < hit_distance ) {
			hit_distance = this_distance;
			result = w->materials[ current_plane.material_index ].color;
		}
	}

	for( u32 sphere_index=0; sphere_index < w->sphere_count; sphere_index++ ) {
		sphere current_sphere = w->spheres[sphere_index];
		
		f32 this_distance = ray_intersects_sphere( ray_origin, ray_direction, current_sphere );
		if( this_distance > 0 && this_distance < hit_distance ) {
			hit_distance = this_distance;
			result = w->materials[ current_sphere.material_index ].color;
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

int main() {

	image_rgba img = allocate_image( 800, 600 );
	
	material mat[4];
	mat[0].color = V3(0.1f, 0.1f, 0.1f);
	mat[1].color = V3(1,0,0);
	mat[2].color = V3(0,0,1);
	mat[3].color = V3(0,0.5f,0.5f);
	
	plane p;
	p.N = V3( 0, 0, 1.0f );
	p.d = 0;
	p.material_index = 1;
	
	sphere spheres[1];
	spheres[0].r = 1;
	spheres[0].P = V3(0,0,0);
	spheres[0].material_index = 2;
	// spheres[1].r = 1;
	// spheres[1].P = V3(1.0f,2.0f,2.0f);
	// spheres[1].material_index = 3;
	
	world w;
	w.material_count = 4;
	w.materials = mat;
	w.plane_count = 1;
	w.planes = &p;
	w.sphere_count = 1;
	w.spheres = spheres;
	
	// let's do this!
	// so we pass a ray from every point on the image in the camera direction onto the world and see what we hit.
	// if we hit nothing, we use the mat[0] to set the pixel
	
	v3 camera_p = V3(0, -10, 1);
	// look at the origin
	v3 camera_z = NOZ(camera_p);
	v3 camera_x = NOZ(cross( V3(0,0,1), camera_z ));
	v3 camera_y = NOZ(cross( camera_z, camera_x ));
	
	f32 film_dist = 1.0f;
	f32 film_w = 1.0f;
	f32 film_h = 1.0f;
	f32 half_film_w = 0.5f * film_w;
	f32 half_film_h = 0.5f * film_h;
	v3 film_center = camera_p - film_dist * camera_z;
	
	u32* out = img.pixels;
	// convert these from image space to film space which is a square from -1 to 1
	for( u32 y=0; y<img.height; y++ ) {
		if( y % 8 == 0 ) {
			printf("\rRaycasting %.2f%%", 100.0 * (f32)y / (f32)img.height );
		}
		f32 film_y = -1.0f + 2.0f * ( (f32)y / (f32)img.height );
		for( u32 x=0; x<img.width; x++ ) {
			
			
			f32 film_x = -1.0f + 2.0f * ( (f32)x / (f32)img.width );
			
			v3 film_p = film_center + film_x * half_film_w * camera_x + film_y * half_film_h * camera_y;
			
			v3 ray_origin = camera_p;
			v3 ray_direction = NOZ( film_p - camera_p );
			
			v3 color = ray_cast( &w, ray_origin, ray_direction );
			u32 color_bmp = rgbapack_4x8( color );
			*out++ = color_bmp;
		}
	}
	
	
	write_bmp( "test.bmp", img );
	
	return 0;
}

