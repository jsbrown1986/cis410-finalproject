/*=========================================================================

			Final Project for CIS 410 (W18)
			GetBin, TriInterp, and Main code 
			implemented by Jacob Brown 3/22/2018

			Produces an image of 1024 x 1024 with
			1024 samples per ray

===========================================================================*/
#include <iostream>
#include <cmath>

#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkSmartPointer.h>
#include <vtkRectilinearGridReader.h>

#define DIM_MAX 1024
#define SAMPLE_MAX 1024
#define PI 3.14159265

using namespace std;

struct Camera {
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
};

struct TransferFunction
{
    double          min;
    double          max;
    int             numBins;
    unsigned char  *colors;  // size is 3*numBins
    double         *opacities; // size is numBins

    void ApplyTransferFunction(double value, unsigned char *RGB, double &opacity)
    {
        int bin = GetBin(value);

		if (bin != -1) {
			RGB[0] = colors[3 * bin + 0];
			RGB[1] = colors[3 * bin + 1];
			RGB[2] = colors[3 * bin + 2];
			opacity = opacities[bin];
		} else {
			RGB[0] = 0;
			RGB[1] = 0;
			RGB[2] = 0;
			opacity = 0;
		}
    }

	// Determines which bin to return
	int GetBin(double value) {

		// Checks if value's within range, bails if not
		if (value < min || value > max)
			return -1;

		// "Quickly" compares value to bins to see
		// which one to return
		double adj = (max - min) / (double)numBins;

		for (int i = 0; i <= numBins; i++) {
			double low = min + (i * adj);
			double high = min + ((i + 1) * adj);

			if (value >= low && value < high)
				return i;
		}

		// Returns -1 if something goes wrong
		return -1;
	}
};

TransferFunction SetupTransferFunction(void){
	int  i;

	TransferFunction rv;
	rv.min = 10;
	rv.max = 15;
	rv.numBins = 256;
	rv.colors = new unsigned char[3 * 256];
	rv.opacities = new double[256];
	unsigned char charOpacity[256] = {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 14, 14, 14, 14, 14, 14, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 5, 4, 3, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17, 17, 17, 17, 17, 17, 16, 16, 15, 14, 13, 12, 11, 9, 8, 7, 6, 5, 5, 4, 3, 3, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 16, 18, 20, 22, 24, 27, 29, 32, 35, 38, 41, 44, 47, 50, 52, 55, 58, 60, 62, 64, 66, 67, 68, 69, 70, 70, 70, 69, 68, 67, 66, 64, 62, 60, 58, 55, 52, 50, 47, 44, 41, 38, 35, 32, 29, 27, 24, 22, 20, 20, 23, 28, 33, 38, 45, 51, 59, 67, 76, 85, 95, 105, 116, 127, 138, 149, 160, 170, 180, 189, 198, 205, 212, 217, 221, 223, 224, 224, 222, 219, 214, 208, 201, 193, 184, 174, 164, 153, 142, 131, 120, 109, 99, 89, 79, 70, 62, 54, 47, 40, 35, 30, 25, 21, 17, 14, 12, 10, 8, 6, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
	};

	for (i = 0; i < 256; i++)
		rv.opacities[i] = charOpacity[i] / 255.0;
	const int numControlPoints = 8;
	unsigned char controlPointColors[numControlPoints * 3] = {
		71, 71, 219, 0, 0, 91, 0, 255, 255, 0, 127, 0,
		255, 255, 0, 255, 96, 0, 107, 0, 0, 224, 76, 76
	};
	double controlPointPositions[numControlPoints] = { 0, 0.143, 0.285, 0.429, 0.571, 0.714, 0.857, 1.0 };
	for (i = 0; i < numControlPoints - 1; i++)
	{
		int start = controlPointPositions[i] * rv.numBins;
		int end = controlPointPositions[i + 1] * rv.numBins + 1;
		//cerr << "Working on " << i << "/" << i + 1 << ", with range " << start << "/" << end << endl;
		if (end >= rv.numBins)
			end = rv.numBins - 1;
		for (int j = start; j <= end; j++)
		{
			double proportion = (j / (rv.numBins - 1.0) - controlPointPositions[i]) / (controlPointPositions[i + 1] - controlPointPositions[i]);
			if (proportion < 0 || proportion > 1.)
				continue;
			for (int k = 0; k < 3; k++)
				rv.colors[3 * j + k] = proportion*(controlPointColors[3 * (i + 1) + k] - controlPointColors[3 * i + k])
				+ controlPointColors[3 * i + k];
		}
	}

	return rv;
}

Camera SetupCamera(void){
	Camera rv;
	rv.focus[0] = 0;
	rv.focus[1] = 0;
	rv.focus[2] = 0;
	rv.up[0] = 0;
	rv.up[1] = -1;
	rv.up[2] = 0;
	rv.angle = 30;
	rv.near = 7.5e+7;
	rv.far = 1.4e+8;
	rv.position[0] = -8.25e+7;
	rv.position[1] = -3.45e+7;
	rv.position[2] = 3.35e+7;

	return rv;
}

// WriteImage - handy function for creating a .png file
void WriteImage(vtkImageData *img, const char *filename) {
	std::string full_filename = filename;
	full_filename += ".png";
	vtkPNGWriter *writer = vtkPNGWriter::New();
	writer->SetInputData(img);
	writer->SetFileName(full_filename.c_str());
	writer->Write();
	writer->Delete();
}

// Trilinear Interpolation - based on Wikipedia's article on the subject
float TriInterp(float x, float x0, float x1, float y, float y0, float y1, float z, float z0, float z1, 
	float c000, float c100, float c001, float c101, float c010, float c110, float c011, float c111) {

	float xd = (x - x0) / (x1 - x0);
	float yd = (y - y0) / (y1 - y0);
	float zd = (z - z0) / (z1 - z0);
	
	float c00 = c000 * (1.0f - xd) + (c100 * xd);
	float c01 = c001 * (1.0f - xd) + (c101 * xd);
	float c10 = c010 * (1.0f - xd) + (c110 * xd);
	float c11 = c011 * (1.0f - xd) + (c111 * xd);

	float c0 = c00 * (1.0f - yd) + (c10 * yd);
	float c1 = c01 * (1.0f - yd) + (c11 * yd);

	return c0 * (1.0f - zd) + (c1 * zd);
}

// GetPointIndex - returns the index in a 3D space
int GetPointIndex(const int *idx, const int *dims) {
	return idx[2] * dims[0] * dims[1] + idx[1] * dims[0] + idx[0];
}

int main() {
	TransferFunction tf = SetupTransferFunction();

	// Camera
	Camera cam = SetupCamera();

	// Reads dataset
	vtkSmartPointer<vtkRectilinearGridReader> rdr =
		vtkSmartPointer<vtkRectilinearGridReader>::New();
	rdr->SetFileName("astro512_ascii.vtk");
	rdr->Update();

	int dims[3];
	vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *)rdr->GetOutput();
	rgrid->GetDimensions(dims);

	float *X = (float *)rgrid->GetXCoordinates()->GetVoidPointer(0);
	float *Y = (float *)rgrid->GetYCoordinates()->GetVoidPointer(0);
	float *Z = (float *)rgrid->GetZCoordinates()->GetVoidPointer(0);

	// For use with "astro512_ascii.vtk"
	float *F = (float *)rgrid->GetPointData()->GetArray("node_sMD")->GetVoidPointer(0);

	// Initializes image to be produced
	vtkImageData * image = vtkImageData::New();
	image->SetDimensions(DIM_MAX, DIM_MAX, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

	// Initializes buffer
	unsigned char * buffer;
	buffer = (unsigned char *)image->GetScalarPointer(0, 0, 0);

	for (int w = 0; w < (3 * DIM_MAX * DIM_MAX); w++)
		buffer[w] = 0;

	// Get "look" (focus - position)
	double look[3];
	look[0] = cam.focus[0] - cam.position[0];
	look[1] = cam.focus[1] - cam.position[1];
	look[2] = cam.focus[2] - cam.position[2];

	// Finds the cross-product of "look" and "up
	double look_up_cross[3];
	look_up_cross[0] = (look[1] * cam.up[2]) - (look[2] * cam.up[1]);
	look_up_cross[1] = (look[2] * cam.up[0]) - (look[0] * cam.up[2]);
	look_up_cross[2] = (look[0] * cam.up[1]) - (look[1] * cam.up[0]);

	// Normalizing look_up_cross
	double mag_look_up_cross = sqrt((look_up_cross[0] * look_up_cross[0]) + (look_up_cross[1] * look_up_cross[1]) + (look_up_cross[2] * look_up_cross[2]));
	
	double ru[3];
	ru[0] = look_up_cross[0] / mag_look_up_cross;
	ru[1] = look_up_cross[1] / mag_look_up_cross;
	ru[2] = look_up_cross[2] / mag_look_up_cross;

	// "Look... are you cross??"
	double look_ru_cross[3];
	look_ru_cross[0] = (look[1] * ru[2]) - (look[2] * ru[1]);
	look_ru_cross[1] = (look[2] * ru[0]) - (look[0] * ru[2]);
	look_ru_cross[2] = (look[0] * ru[1]) - (look[1] * ru[0]);

	// Normalizing look_ru_cross
	double mag_look_ru_cross = sqrt((look_ru_cross[0] * look_ru_cross[0]) + (look_ru_cross[1] * look_ru_cross[1]) + (look_ru_cross[2] * look_ru_cross[2]));
	
	double rv[3];
	rv[0] = look_ru_cross[0] / mag_look_ru_cross;
	rv[1] = look_ru_cross[1] / mag_look_ru_cross;
	rv[2] = look_ru_cross[2] / mag_look_ru_cross;

	// Convert radians to degrees with tan function
	double rad_to_degrees = tan((cam.angle / 2.0f) * (PI / 180.0f));

	// Delta RX
	double delta_rx[3];
	delta_rx[0] = ((2.0f * rad_to_degrees) / DIM_MAX) * ru[0];
	delta_rx[1] = ((2.0f * rad_to_degrees) / DIM_MAX) * ru[1];
	delta_rx[2] = ((2.0f * rad_to_degrees) / DIM_MAX) * ru[2];

	// Delta RY
	double delta_ry[3];
	delta_ry[0] = ((2.0f * rad_to_degrees) / DIM_MAX) * rv[0];
	delta_ry[1] = ((2.0f * rad_to_degrees) / DIM_MAX) * rv[1];
	delta_ry[2] = ((2.0f * rad_to_degrees) / DIM_MAX) * rv[2];

	// Find magnitude of look
	double mag_look = sqrt((look[0] * look[0]) + (look[1] * look[1]) + (look[2] * look[2]));

	// move me later
	float sample_map[SAMPLE_MAX];

	// Loops through each pixel, determines ray directions,
	// steps through rays, and colors each pixel
	for (int i = 0; i < DIM_MAX; i++){
		for (int j = 0; j < DIM_MAX; j++){

			// Ray direction
			double rd[3];
			rd[0] = (look[0] / mag_look) + ((((2.0f * i) + 1.0f - DIM_MAX) / 2.0f) * delta_rx[0]) + ((((2.0f * j) + 1.0f - DIM_MAX) / 2.0f) * delta_ry[0]);
			rd[1] = (look[1] / mag_look) + ((((2.0f * i) + 1.0f - DIM_MAX) / 2.0f) * delta_rx[1]) + ((((2.0f * j) + 1.0f - DIM_MAX) / 2.0f) * delta_ry[1]);
			rd[2] = (look[2] / mag_look) + ((((2.0f * i) + 1.0f - DIM_MAX) / 2.0f) * delta_rx[2]) + ((((2.0f * j) + 1.0f - DIM_MAX) / 2.0f) * delta_ry[2]);

			// Ray Sampling
			double near_adj[3];
			near_adj[0] = cam.near * rd[0];
			near_adj[1] = cam.near * rd[1];
			near_adj[2] = cam.near * rd[2];

			double far_adj[3];
			far_adj[0] = cam.far * rd[0];
			far_adj[1] = cam.far * rd[1];
			far_adj[2] = cam.far * rd[2];

			double cur_pos[3];
			cur_pos[0] = cam.position[0] + near_adj[0];
			cur_pos[1] = cam.position[1] + near_adj[1];
			cur_pos[2] = cam.position[2] + near_adj[2];

			double delta[3];
			delta[0] = (far_adj[0] - near_adj[0]) / (SAMPLE_MAX - 1);
			delta[1] = (far_adj[1] - near_adj[1]) / (SAMPLE_MAX - 1);
			delta[2] = (far_adj[2] - near_adj[2]) / (SAMPLE_MAX - 1);

			int idx[3];

			// Variables for color calculations
			unsigned char RGB[3];
			RGB[0] = 0;
			RGB[1] = 0;
			RGB[2] = 0;

			double fr = 0.0f;
			double fg = 0.0f;
			double fb = 0.0f;
			double fa = 0.0f;

			double br = 0.0f;
			double bg = 0.0f;
			double bb = 0.0f;
			double ba = 0.0f;

			double opacity = 0.0f;

			double r = 0.0f;
			double g = 0.0f;
			double b = 0.0f;
			double a = 0.0f;

			float sample;

			// Runs through each step, determines samples, then
			// colors each pixel
			for (int p = 0; p < SAMPLE_MAX; p++) {

				// Checks to see if it's within X, Y, Z range
				if (cur_pos[0] >= X[0] && cur_pos[0] < X[dims[0] - 1] &&
					cur_pos[1] >= Y[0] && cur_pos[1] < Y[dims[1] - 1] &&
					cur_pos[2] >= Z[0] && cur_pos[2] < Z[dims[2] - 1]) {

					// Check X
					for (int a = 0; a < dims[0] - 1; a++) {
						if (cur_pos[0] >= X[a] && cur_pos[0] < X[a + 1]) {
							idx[0] = a;
							break;
						}
					}

					// Check Y
					for (int b = 0; b < dims[1] - 1; b++) {
						if (cur_pos[1] >= Y[b] && cur_pos[1] < Y[b + 1]) {
							idx[1] = b;
							break;
						}
					}

					// Check Z
					for (int c = 0; c < dims[2] - 1; c++) {
						if (cur_pos[2] >= Z[c] && cur_pos[2] < Z[c + 1]) {
							idx[2] = c;
							break;
						}
					}

					// Grabs logical cell indices
					int c000_log[3], c100_log[3], c010_log[3], c110_log[3], c001_log[3], c101_log[3], c011_log[3], c111_log[3];
					float c000_f, c100_f, c010_f, c110_f, c001_f, c101_f, c011_f, c111_f;
					
					c000_log[0] = idx[0];
					c000_log[1] = idx[1];
					c000_log[2] = idx[2];

					c100_log[0] = idx[0] + 1;
					c100_log[1] = idx[1];
					c100_log[2] = idx[2];

					c010_log[0] = idx[0];
					c010_log[1] = idx[1] + 1;
					c010_log[2] = idx[2];

					c110_log[0] = idx[0] + 1;
					c110_log[1] = idx[1] + 1;
					c110_log[2] = idx[2];

					c001_log[0] = idx[0];
					c001_log[1] = idx[1];
					c001_log[2] = idx[2] + 1;

					c101_log[0] = idx[0] + 1;
					c101_log[1] = idx[1];
					c101_log[2] = idx[2] + 1;

					c011_log[0] = idx[0];
					c011_log[1] = idx[1] + 1;
					c011_log[2] = idx[2] + 1;

					c111_log[0] = idx[0] + 1;
					c111_log[1] = idx[1] + 1;
					c111_log[2] = idx[2] + 1;

					// Grabs field values at each index
					c000_f = F[GetPointIndex(c000_log, dims)];
					c100_f = F[GetPointIndex(c100_log, dims)];
					c010_f = F[GetPointIndex(c010_log, dims)];
					c110_f = F[GetPointIndex(c110_log, dims)];

					c001_f = F[GetPointIndex(c001_log, dims)];
					c101_f = F[GetPointIndex(c101_log, dims)];
					c011_f = F[GetPointIndex(c011_log, dims)];
					c111_f = F[GetPointIndex(c111_log, dims)];

					// Sample at this location
					sample = TriInterp(cur_pos[0], X[idx[0]], X[idx[0] + 1], cur_pos[1], Y[idx[1]], Y[idx[1] + 1], cur_pos[2], Z[idx[2]], Z[idx[2] + 1],
							c000_f, c100_f, c001_f, c101_f, c010_f, c110_f, c011_f, c111_f);
				} else {
					sample = 0.0f;
				}

				// Increments position on ray
				cur_pos[0] += delta[0];
				cur_pos[1] += delta[1];
				cur_pos[2] += delta[2];

				// Begin color stuff
				fr = r;
				fg = g;
				fb = b;
				fa = a;

				// Transfer function and opacity correction
				tf.ApplyTransferFunction(sample, RGB, opacity);
				opacity = 1.0f - pow((1.0f - opacity), (500.0f / SAMPLE_MAX));

				br = (double)RGB[0] / 255.0f;
				bg = (double)RGB[1] / 255.0f;
				bb = (double)RGB[2] / 255.0f;
				ba = opacity;

				(br > 0.0f) ? r = fr + (1.0f - fa) * ba * br : r = fr;
				(bg > 0.0f) ? g = fg + (1.0f - fa) * ba * bg : g = fg;
				(bb > 0.0f) ? b = fb + (1.0f - fa) * ba * bb : b = fb;
				(ba > 0.0f) ? a = fa + (1.0f - fa) * ba : a = fa;

				// If the opacity is solid, break
				if (a == 1.0f)
					break;
			}

			// Sends color info to buffer
			int offset = 3 * (j * DIM_MAX + i);
			buffer[offset + 0] = (int)(r * 255);
			buffer[offset + 1] = (int)(g * 255);
			buffer[offset + 2] = (int)(b * 255);
		}
	}

	WriteImage(image, "jbrown_final");
	return 0;
}