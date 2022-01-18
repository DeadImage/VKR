#define _CRT_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <string>
#include "wavelets.h"
#include <iostream>

#define CHANNEL_NUM 1

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace std;

int main() {

	int width, height, bpp;
	uint8_t* imageSam = stbi_load("imageSam.jpg", &width, &height, &bpp, CHANNEL_NUM);
	int width1, height1, bpp1;
	uint8_t* imageLena = stbi_load("imageLena.jpg", &width1, &height1, &bpp1, CHANNEL_NUM);

	uint8_t* copyOfSam = new uint8_t [width*height];
	uint8_t* copyOfLena = new uint8_t[width1 * height1];
	for (int i = 0; i < width * height; i++)
		copyOfSam[i] = imageSam[i];
	for (int i = 0; i < width1 * height1; i++)
		copyOfLena[i] = imageLena[i];


	dwt_2d_lifting(copyOfSam, height, width, 1, "haar");
	dwt_2d_lifting(copyOfLena, height1, width1, 1, "haar");

	stbi_write_jpg("Lena_lifting_Haar.jpg", width1, height1, CHANNEL_NUM, copyOfLena, 100);
	stbi_write_jpg("Sam_lifting_Haar.jpg", width, height, CHANNEL_NUM, copyOfSam, 100);

	for (int i = 0; i < width * height; i++)
		copyOfSam[i] = imageSam[i];
	for (int i = 0; i < width1 * height1; i++)
		copyOfLena[i] = imageLena[i];

	dwt_2d_lifting(copyOfSam, height, width, 1, "db2");
	dwt_2d_lifting(copyOfLena, height1, width1, 1, "db2");

	stbi_write_jpg("Lena_lifting_db2.jpg", width1, height1, CHANNEL_NUM, copyOfLena, 100);
	stbi_write_jpg("Sam_lifting_db2.jpg", width, height, CHANNEL_NUM, copyOfSam, 100);

	for (int i = 0; i < width * height; i++)
		copyOfSam[i] = imageSam[i];
	for (int i = 0; i < width1 * height1; i++)
		copyOfLena[i] = imageLena[i];

	dwt_haar_2d(copyOfSam, height, width, 1);
	dwt_haar_2d(copyOfLena, height1, width1, 1);

	stbi_write_jpg("Lena_filtering_Haar.jpg", width1, height1, CHANNEL_NUM, copyOfLena, 100);
	stbi_write_jpg("Sam_filtering_Haar.jpg", width, height, CHANNEL_NUM, copyOfSam, 100);

	return 0;
}