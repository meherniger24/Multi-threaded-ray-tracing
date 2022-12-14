#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<cmath>
#include<thread>
#include<chrono>
#include<ctime>

using namespace std;


struct vec3 {
	float x, y, z;
	float len;
	vec3(float X, float Y, float Z)
		: x(X), y(Y), z(Z)
	{
		len = sqrt(X * X + Y * Y + Z * Z);
	}
	vec3()
		: x(0), y(0), z(0), len(0)
	{}
	vec3 operator-(vec3 tar) {
		return vec3(x - tar.x, y - tar.y, z - tar.z);
	}
	vec3 operator+(vec3 tar) {
		return vec3(x + tar.x, y + tar.y, z + tar.z);
	}
	// dot (inner) product of two vectors
	float dot(vec3& tar) {
		return (x * tar.x) + (y * tar.y) + (z * tar.z);
	}

	vec3 operator*(vec3 tar) {
		return vec3(x * tar.x, y * tar.y, z * tar.z);
	}
	vec3 operator*(float t) {
		return vec3(x * t, y * t, z * t);
	}
	vec3 operator/(float t) {
		return vec3(x / t, y / t, z / t);
	}
	vec3 normalize() {
		return vec3(x / len, y / len, z / len);
	}
	
};



struct color {
	float r, g, b;
	color(float r, float g, float b)
		: r(r), g(g), b(b)
	{}
	color()
		: r(0), g(0), b(0)
	{}
	color operator*(float x) {
		r *= x;
		g *= x;
		b *= x;
		return *this;
	}
	color operator*(color clr) {
		r *= clr.r;
		g *= clr.g;
		b*= clr.b;
		return *this;
	}
	color operator+(color clr) {
		color result;
		result.r = r * clr.r;
		result.b = b * clr.b;
		result.g = g * clr.g;
		return result;
	}
};


	rays(vec3 p, vec3 vec)
		: view_point(p), v(vec)
	{}
};



struct light {
	vec3 center;
	color clr;
};



//when the input ray intersects the  sphere, the ouput color is an array of color r, g, b

color intersection(rays ray, vector<sphere> s, vector<light> lightt) {
	color rgb(1.0f, 1.0f, 1.0f);
	//color light_color;
	float t = 0.0f;
	float discriminant;
	float a, b, c;
	float minimum = (float)INT_MAX;							


	for (int i = 0; i < s.size(); i++) {
		sphere sph = s.at(i);
		vec3 distance = ray.view_point - sph.center;
		a = ray.v.dot(ray.v);
		b = ray.v.dot(distance) * 2;
		c = distance.dot(distance) - (sph.radius * sph.radius);
		discriminant = (b * b) - (4 * a * c);
		t = (-b - sqrt(discriminant)) / (2 * a);
		if (discriminant >= 0 && t < minimum) {
			minimum = t;
			float intensity = 0.0f;			
			rays light_ray;
			light_ray.view_point = ray.v * t;

			//calculate unit vector from view_point to center of light
			for (int i = 0; i < lightt.size(); i++) {
				light lit = lightt.at(i);
				light_ray.v = (lit.center - sph.center).normalize();

				//calculating illumination
				intensity += (light_ray.v.dot(Normal) > 0) ? light_ray.v.dot(Normal) : 0;
			}
			intensity /= lightt.size();
			rgb = sph.clr * intensity;
		}
	}
	return rgb;


};


//reading_sphere
vector<sphere> RSphere(string s) {
	ifstream readfile(s);
	vector<sphere> spheres;
	float x, y, z;
	while (readfile) {
		sphere sph;
		readfile >> x >> y >> z;
		sph.center.x = x;
		sph.center.y = y;
		sph.center.z = z;
		readfile >> sph.radius;
		readfile >> sph.clr.b >> sph.clr.g >> sph.clr.r;
		spheres.push_back(sph);
	}
	return spheres;
}

//reading_light
vector<light> RLight(string s) {
	ifstream rfile(s);
	vector<light> Lights;
	float x, y, z;
	while (rfile) {
		light light1;
		rfile >> x >> y >> z;
		light1.center.x = x;
		light1.center.y = y;
		light1.center.z = z;
		rfile >> light1.clr.r >> light1.clr.g >> light1.clr.b;
		Lights.push_back(light1);
	}
	return Lights;
	
}


// function_for saving targa image

void write_image(string filename, char* bytes, short N) {
	std::ofstream outfile;
	outfile.open(filename, std::ios::binary | std::ios::out);	// open a binary file
	outfile.put(0);						// id length (field 1)
	outfile.put(0);						// color map type (field 2)
	outfile.put(2);						// image_type (field 3)
	outfile.put(0); outfile.put(0);		// color map field entry index (field 4)
	outfile.put(0); outfile.put(0);		// color map length (field 4)
	outfile.put(0);				// color map entry size (field 4)
	outfile.put(0); outfile.put(0);		// x origin (field 5)
	outfile.put(0); outfile.put(0);		// y origin (field 5)
	outfile.write((char*)&N, 2);		// image width (field 5)
	outfile.write((char*)&N, 2);		// image height (field 5)
	outfile.put(24);				// pixel depth (field 5)
	outfile.put(0);				// image descriptor (field 5)
	outfile.write(bytes, N * N * 3);		// write the image data
	outfile.close();				// close the file
}


// raytraced_image creating

void raytracedImage(char* image, vector<sphere> spheres, vector<light> lights, short N, int start, int end) {
	float z = 1.5;
	float dx = 2.0f / (N - 1);			//change scale
	float dy = 2.0f / (N - 1);
	float sx = -1.0f;		            //starting point
	float sy = -1.0f;
	for (short i = start; i < end; i++) {
		for (short j = 0; j < N; j++) {
			rays ray;
			ray.v = direction;
			ray.view_point = view;
			color pixel = intersection(ray, spheres, lights);
			image[i * N * 3 + j * 3 + 0] = pixel.r * 255;
			image[i * N * 3 + j * 3 + 1] = pixel.g * 255;
			image[i * N * 3 + j * 3 + 2] = pixel.b * 255;
		}
	}

}


int main(int argc, char** argv[])
{
	
	short N = 1024;
	
	char* image = (char*)malloc(N * N * 3);
	vector<sphere> spheres;
	vector<light> lights;



	spheres = RSphere("spheres.txt");
	lights = RLight("lights.txt");

	// number of threads 
	int numThreads = 100; 


	vector<std::thread> threads;

	// profiling
	time_t T1 = time(NULL);
	
	for (int i = 0; i < numThreads; i++) {
		int start = (double)i / (double)numThreads * N;
		int stop = (double)(i + 1) / (double)numThreads * N;
		threads.push_back(std::thread(raytracedImage, image, spheres, lights, N, start, stop));
	}

	for (auto& thread : threads) {
		thread.join();
	}

	std::time_t T2 = time(NULL);
	cout << "Time to compute: " << (T2 - T1) << endl;
	

	write_image("spheres.tga", image, N);
	return 0;
}
