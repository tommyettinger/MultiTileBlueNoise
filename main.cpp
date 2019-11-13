#define _CRT_SECURE_NO_WARNINGS

#define THRESHOLD_SAMPLES() 11 // the number of samples for threshold testing.

#include <vector>
#include <stdint.h>

//#include "scoped_timer.h"
#include "generatebn_void_cluster.h"


int main(int argc, char** argv)
{
	// generate blue noise using void and cluster
	{
		static size_t c_width = 128;

		std::vector<uint8_t> noise;

		{
//			ScopedTimer timer("Blue noise by void and cluster");

			GenerateBN_Void_Cluster(noise, c_width, false, "out/blue");
		}
		//TestNoise(noise, c_width, "out/blueVC_%i");
	}

    //// generate blue noise using void and cluster but using mitchell's best candidate instead of initial binary pattern and phase 1
    //{
    //    static size_t c_width = 256;

    //    std::vector<uint8_t> noise;

    //    {
    //        ScopedTimer timer("Blue noise by void and cluster with Mitchells best candidate");
    //        GenerateBN_Void_Cluster(noise, c_width, true, "out/blueVC_1M");
    //    }

    //    TestNoise(noise, c_width, "out/blueVC_1M");
    //}

    //// load a blue noise texture
    //{
    //    int width, height, channels;

    //    std::vector<uint8_t> noise;

    //    {
    //        ScopedTimer timer("Blue noise by void and cluster from loading the texture");
    //        uint8_t* image = stbi_load("bluenoise256.png", &width, &height, &channels, 0);

    //        noise.reserve(width*height);
    //        for (int i = 0; i < width*height; ++i)
    //            noise.push_back(image[i*channels]);

    //        stbi_image_free(image);
    //    }

    //    TestNoise(noise, width, "out/blueVC_2");
    //}
	
	//system("pause");

    return 0;
}