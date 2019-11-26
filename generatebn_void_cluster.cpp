#define _CRT_SECURE_NO_WARNINGS

#include "generatebn_void_cluster.h"
#include "whitenoise.h"
#include "convert.h"

#include "dft.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

#include "histogram.h"
#include "image.h"

#include "scoped_timer.h"

static const float c_sigma = 1.9f;// 3.5f;// 3.6180339887498949f;// +0.7548776662466927f;// 1.9f;// 1.5f;
static const float c_2sigmaSquared = 2.0f * c_sigma * c_sigma;
static const int c_3sigmaint = int(ceil(c_sigma * 3.0f));

static void SaveLUTImage(const std::vector<bool>& binaryPattern, std::vector<float>& LUT, size_t width, const char* fileName)
{
    // get the LUT min and max
    float LUTMin = LUT[0];
    float LUTMax = LUT[0];
    for (float f : LUT)
    {
        LUTMin = std::min(LUTMin, f);
        LUTMax = std::max(LUTMax, f);
    }

    size_t c_scale = 4;

    std::vector<uint8_t> image(width*width * c_scale*c_scale * 3);
    for (size_t index = 0; index < width*width*c_scale*c_scale; ++index)
    {
        size_t x = (index % (width * c_scale)) / c_scale;
        size_t y = index / (width * c_scale * c_scale);

        float percent = (LUT[y*width + x] - LUTMin) / (LUTMax - LUTMin);
        uint8_t value = FromFloat<uint8_t>(percent);

        image[index * 3 + 0] = value;
        image[index * 3 + 1] = value;
        image[index * 3 + 2] = value;

        if (binaryPattern[y*width + x])
        {
            image[index * 3 + 0] = 0;
            image[index * 3 + 1] = 255;
            image[index * 3 + 2] = 0;
        }
    }
    stbi_write_png(fileName, int(width*c_scale), int(width*c_scale), 3, image.data(), 0);
}

#if 1

template <bool CLUSTER>
static bool FindWinnerLUT(const std::vector<float>& LUT, const std::vector<bool>& binaryPattern, size_t width, int &bestPixelX, int& bestPixelY, std::mt19937& rng)
{
    float bestValue = CLUSTER ? -FLT_MAX : FLT_MAX;
    std::vector<size_t> bestIndices;
    for (size_t index = 0, count = LUT.size(); index < count; ++index)
    {
        if (binaryPattern[index] == CLUSTER)
        {
            if (LUT[index] == bestValue)
            {
                bestIndices.push_back(index);
            }
            else if ((CLUSTER == true && LUT[index] > bestValue) || (CLUSTER == false && LUT[index] < bestValue))
            {
                bestValue = LUT[index];
                bestIndices.clear();
                bestIndices.push_back(index);
            }
        }
    }

    if (bestIndices.size() == 0)
        return false;

    size_t bestIndex = bestIndices[0];

    // can randomize the winners
    /*
    if (bestIndices.size() > 1)
    {
        std::uniform_int_distribution<size_t> dist(0, bestIndices.size() - 1);
        bestIndex = bestIndices[dist(rng)];
    }
    */

    bestPixelX = int(bestIndex % width);
    bestPixelY = int(bestIndex / width);

    return true;
}

static bool FindTightestClusterLUT(const std::vector<float>& LUT, const std::vector<bool>& binaryPattern, size_t width, int &bestPixelX, int& bestPixelY, std::mt19937& rng)
{
    return FindWinnerLUT<true>(LUT, binaryPattern, width, bestPixelX, bestPixelY, rng);
}

static bool FindLargestVoidLUT(const std::vector<float>& LUT, const std::vector<bool>& binaryPattern, size_t width, int &bestPixelX, int& bestPixelY, std::mt19937& rng)
{
    return FindWinnerLUT<false>(LUT, binaryPattern, width, bestPixelX, bestPixelY, rng);
}

#else
static bool FindTightestClusterLUT(const std::vector<float>& LUT, const std::vector<bool>& binaryPattern, size_t width, int &bestPixelX, int& bestPixelY, std::mt19937& rn)
{
    float bestValue = -FLT_MAX;
    size_t bestIndex = ~size_t(0);
    for (size_t index = 0, count = LUT.size(); index < count; ++index)
    {
        if (binaryPattern[index] && LUT[index] > bestValue)
        {
            bestValue = LUT[index];
            bestIndex = index;
        }
    }

    if (bestIndex == ~size_t(0))
        return false;

    bestPixelX = int(bestIndex % width);
    bestPixelY = int(bestIndex / width);

    return true;
}

static bool FindLargestVoidLUT(const std::vector<float>& LUT, const std::vector<bool>& binaryPattern, size_t width, int &bestPixelX, int& bestPixelY, std::mt19937& rn)
{
    float bestValue = FLT_MAX;
    size_t bestIndex = ~size_t(0);
    for (size_t index = 0, count = LUT.size(); index < count; ++index)
    {
        if (!binaryPattern[index] && LUT[index] < bestValue)
        {
            bestValue = LUT[index];
            bestIndex = index;
        }
    }

    if (bestIndex == ~size_t(0))
        return false;

    bestPixelX = int(bestIndex % width);
    bestPixelY = int(bestIndex / width);

    return true;
}
#endif

static void WriteLUTValue(std::vector<float>& LUT, size_t width, bool value, int basex, int basey)
{
    #pragma omp parallel for
    for (int y = 0; y < width; ++y)
    {
        float disty = abs(float(y) - float(basey));

		if (disty > float(width) * 0.5f)
			disty = float(width) - disty;

//		if (float(basey) < float(width) * 0.125f)
//		{
//			disty = std::min(float(y) + float(basey), disty);
//			disty = std::min(abs(float(width) * 1.0f - y) + float(basey), disty);
//			if(y >= float(width) * 0.125f && y < float(width) * 0.25f)
//				disty = std::min(abs(float(width) * 0.25f - y) + float(basey), disty);
//		}
//		else if (float(basey) >= float(width) * 0.125f && float(basey) < float(width) * 0.25f)
//		{
//			disty = std::min(abs(float(width) * 0.25f - y) + abs(float(width) * 0.25f - basey), disty);
//			if (y < float(width) * 0.125f)
//				disty = std::min(abs(float(width) * 1.0f - y) + abs(float(width) * 0.25f - basey), disty);
//		}
//		else if (float(basey) >= float(width) * 0.25f && float(basey) < float(width) * 0.375f)
//		{
//			disty = std::min(abs(float(width) * 0.25f - y) + abs(float(width) * 0.25f - basey), disty);
//			if (y >= float(width) * 0.875f)
//				disty = std::min(abs(float(width) * 1.0f - y) + abs(float(width) * 0.25f - basey), disty);
//		}
//		else if (float(basey) >= float(width) * 0.375f && float(basey) < float(width) * 0.5f)
//		{
//			disty = std::min(abs(float(width) * 0.5f - y) + abs(float(width) * 0.5f - basey), disty);
//			if (y >= float(width) * 0.75f && y < float(width) * 0.875f)
//				disty = std::min(abs(float(width) * 0.75f - y) + abs(float(width) * 0.5f - basey), disty);
//		}
//		else if (float(basey) >= float(width) * 0.5f && float(basey) < float(width) * 0.625f)
//		{
//			disty = std::min(abs(float(width) * 0.5f - y) + abs(float(width) * 0.5f - basey), disty);
//			if (y >= float(width) * 0.625f && y < float(width) * 0.75f)
//				disty = std::min(abs(float(width) * 0.75f - y) + abs(float(width) * 0.5f - basey), disty);
//		}
//		else if (float(basey) >= float(width) * 0.625f && float(basey) < float(width) * 0.75f)
//		{
//			disty = std::min(abs(float(width) * 0.75f - y) + abs(float(width) * 0.75f - basey), disty);
//			if (y >= float(width) * 0.5f && y < float(width) * 0.625f)
//				disty = std::min(abs(float(width) * 0.5f - y) + abs(float(width) * 0.75f - basey), disty);
//		}
//		else if (float(basey) >= float(width) * 0.75f && float(basey) < float(width) * 0.875f)
//		{
//			disty = std::min(abs(float(width) * 0.75f - y) + abs(float(width) * 0.75f - basey), disty);
//			if (y >= float(width) * 0.375f && y < float(width) * 0.5f)
//				disty = std::min(abs(float(width) * 0.5f - y) + abs(float(width) * 0.75f - basey), disty);
//		}
//		else if (float(basey) >= float(width) * 0.875f)
//		{
//			disty = std::min(abs(float(width) * 1.0f - y) + abs(float(width) * 1.0f - basey), disty);
//			disty = std::min(float(y) + abs(float(width) * 1.0f - basey), disty);
//			if (y >= float(width) * 0.25f && y < float(width) * 0.375f)
//				disty = std::min(abs(float(width) * 0.25f - y) + abs(float(width) * 1.0f - basey), disty);
//		}

        for (size_t x = 0; x < width; ++x)
        {
            float distx = abs(float(x) - float(basex));
			if (distx > float(width) * 0.5f)
				distx = float(width) - distx;

//			if (float(basex) < float(width) * 0.125f)
//			{
//				distx = std::min(float(x) + float(basex), distx);
//				distx = std::min(abs(float(width) * 1.0f - x) + float(basex), distx);
//				if (x >= float(width) * 0.125f && x < float(width) * 0.25f)
//					distx = std::min(abs(float(width) * 0.25f - x) + float(basex), distx);
//			}
//			else if (float(basex) >= float(width) * 0.125f && float(basex) < float(width) * 0.25f)
//			{
//				distx = std::min(abs(float(width) * 0.25f - x) + abs(float(width) * 0.25f - basex), distx);
//				if (x < float(width) * 0.125f)
//					distx = std::min(abs(float(width) * 1.0f - x) + abs(float(width) * 0.25f - basex), distx);
//			}
//			else if (float(basex) >= float(width) * 0.25f && float(basex) < float(width) * 0.375f)
//			{
//				distx = std::min(abs(float(width) * 0.25f - x) + abs(float(width) * 0.25f - basex), distx);
//				if (x >= float(width) * 0.875f)
//					distx = std::min(abs(float(width) * 1.0f - x) + abs(float(width) * 0.25f - basex), distx);
//			}
//			else if (float(basex) >= float(width) * 0.375f && float(basex) < float(width) * 0.5f)
//			{
//				distx = std::min(abs(float(width) * 0.5f - x) + abs(float(width) * 0.5f - basex), distx);
//				if (x >= float(width) * 0.75f && x < float(width) * 0.875f)
//					distx = std::min(abs(float(width) * 0.75f - x) + abs(float(width) * 0.5f - basex), distx);
//			}
//			else if (float(basex) >= float(width) * 0.5f && float(basex) < float(width) * 0.625f)
//			{
//				distx = std::min(abs(float(width) * 0.5f - x) + abs(float(width) * 0.5f - basex), distx);
//				if (x >= float(width) * 0.625f && x < float(width) * 0.75f)
//					distx = std::min(abs(float(width) * 0.75f - x) + abs(float(width) * 0.5f - basex), distx);
//			}
//			else if (float(basex) >= float(width) * 0.625f && float(basex) < float(width) * 0.75f)
//			{
//				distx = std::min(abs(float(width) * 0.75f - x) + abs(float(width) * 0.75f - basex), distx);
//				if (x >= float(width) * 0.5f && x < float(width) * 0.625f)
//					distx = std::min(abs(float(width) * 0.5f - x) + abs(float(width) * 0.75f - basex), distx);
//			}
//			else if (float(basex) >= float(width) * 0.75f && float(basex) < float(width) * 0.875f)
//			{
//				distx = std::min(abs(float(width) * 0.75f - x) + abs(float(width) * 0.75f - basex), distx);
//				if (x >= float(width) * 0.375f && x < float(width) * 0.5f)
//					distx = std::min(abs(float(width) * 0.5f - x) + abs(float(width) * 0.75f - basex), distx);
//			}
//			else if (float(basex) >= float(width) * 0.875f)
//			{
//				distx = std::min(abs(float(width) * 1.0f - x) + abs(float(width) * 1.0f - basex), distx);
//				distx = std::min(float(x) + abs(float(width) * 1.0f - basex), distx);
//				if (x >= float(width) * 0.25f && x < float(width) * 0.375f)
//					distx = std::min(abs(float(width) * 0.25f - x) + abs(float(width) * 1.0f - basex), distx);
//			}


            float distanceSquared = float(distx*distx) + float(disty*disty);
            float energy = exp(-distanceSquared / c_2sigmaSquared) * (value ? 1.0f : -1.0f);
			LUT[y * width + x] += energy;
		}
    }
}

static void MakeLUT(const std::vector<bool>& binaryPattern, std::vector<float>& LUT, size_t width, bool writeOnes)
{
    LUT.clear();
    LUT.resize(width*width, 0.0f);
    for (size_t index = 0; index < width*width; ++index)
    {
        if (binaryPattern[index] == writeOnes)
        {
            int x = int(index % width);
            int y = int(index / width);
            WriteLUTValue(LUT, width, writeOnes, x, y);
        }
    }
}

#if SAVE_VOIDCLUSTER_INITIALBP()

static void SaveBinaryPattern(const std::vector<bool>& binaryPattern, size_t width, const char* baseFileName, int iterationCount, int tightestClusterX, int tightestClusterY, int largestVoidX, int largestVoidY)
{
    size_t c_scale = 4;

    std::vector<uint8_t> binaryPatternImage(width*width * c_scale*c_scale * 3);
    for (size_t index = 0; index < width*width*c_scale*c_scale; ++index)
    {
        size_t x = (index % (width * c_scale)) / c_scale;
        size_t y = index / (width * c_scale * c_scale);

        bool isCluster = (x == tightestClusterX && y == tightestClusterY);
        bool isVoid = (x == largestVoidX && y == largestVoidY);

        if (isCluster == isVoid)
        {
            if (isCluster)
            {
                binaryPatternImage[index * 3 + 0] = 255;
                binaryPatternImage[index * 3 + 1] = 255;
                binaryPatternImage[index * 3 + 2] = 0;
            }
            else
            {
                binaryPatternImage[index * 3 + 0] = binaryPattern[y*width+x] ? 255 : 0;
                binaryPatternImage[index * 3 + 1] = binaryPattern[y*width + x] ? 255 : 0;
                binaryPatternImage[index * 3 + 2] = binaryPattern[y*width + x] ? 255 : 0;
            }
        }
        else if (isCluster)
        {
            binaryPatternImage[index * 3 + 0] = 255;
            binaryPatternImage[index * 3 + 1] = 0;
            binaryPatternImage[index * 3 + 2] = 0;
        }
        else if (isVoid)
        {
            binaryPatternImage[index * 3 + 0] = 0;
            binaryPatternImage[index * 3 + 1] = 255;
            binaryPatternImage[index * 3 + 2] = 0;
        }
    }

    char fileName[256];
    sprintf(fileName, "%s_base_%i.png", baseFileName, iterationCount);
    stbi_write_png(fileName, int(width*c_scale), int(width*c_scale), 3, binaryPatternImage.data(), 0);
}

#endif

struct Point
{
	size_t x;
	size_t y;
};
typedef std::vector<Point> TPoints;
typedef std::vector<TPoints> TPointGrid;

static bool DistanceSqToClosestPoint(const TPoints& points, const Point& point, float& minDistSq, size_t width)
{
	if (points.size() == 0)
		return false;

	// calculate the closest distance from this point to an existing sample
	for (const Point& p : points)
	{
		float distx = std::abs(float(p.x) - float(point.x));
		float disty = std::abs(float(p.y) - float(point.y));

		if (distx > float(width) / 2.0f)
			distx = float(width) - distx;

		if (disty > float(width) / 2.0f)
			disty = float(width) - disty;

		float distSq = distx * distx + disty * disty;
		if (distSq < minDistSq)
			minDistSq = distSq;
	}
	return true;
}

static float DistanceSqToClosestPoint(const TPointGrid& grid, size_t cellCount, size_t cellSize, const Point& point, size_t width)
{
	const int basex = int(point.x / cellSize);
	const int basey = int(point.y / cellSize);

	const int maxRadius = int(cellCount / 2);

	float minDistSq = FLT_MAX;
	bool foundAPoint = false;
	bool didAnExtraRing = false;

	for (int radius = 0; radius <= maxRadius; ++radius)
	{
		// top and bottom rows
		{
			for (int offsetX = -radius; offsetX <= radius; ++offsetX)
			{
				int x = int(basex + offsetX + cellCount) % int(cellCount);

				int offsetY = -radius;
				int y = int(basey + offsetY + cellCount) % int(cellCount);
				foundAPoint |= DistanceSqToClosestPoint(grid[y * cellCount + x], point, minDistSq, width);

				offsetY = radius;
				y = int(basey + offsetY + cellCount) % int(cellCount);
				foundAPoint |= DistanceSqToClosestPoint(grid[y * cellCount + x], point, minDistSq, width);
			}
		}

		// left and right
		{
			for (int offsetY = -radius + 1; offsetY <= radius - 1; ++offsetY)
			{
				int y = int(basey + offsetY + cellCount) % int(cellCount);

				int offsetX = -radius;
				int x = int(basex + offsetX + cellCount) % int(cellCount);
				foundAPoint |= DistanceSqToClosestPoint(grid[y * cellCount + x], point, minDistSq, width);

				offsetX = +radius;
				x = int(basex + offsetX + cellCount) % int(cellCount);
				foundAPoint |= DistanceSqToClosestPoint(grid[y * cellCount + x], point, minDistSq, width);
			}
		}

		// we stop when we've found a point, then do another ring to make sure there isn't something closer to what we found.
		if (foundAPoint)
		{
			if (didAnExtraRing)
				break;
			else
				didAnExtraRing = true;
		}
	}

	return minDistSq;
}

static void AddPointToPointGrid(TPointGrid& grid, size_t cellCount, size_t cellSize, const Point& point)
{
	Point cell;
	cell.x = point.x / cellSize;
	cell.y = point.y / cellSize;
	grid[cell.y * cellCount + cell.x].push_back(point);
}

// This replaces "Initial Binary Pattern" and "Phase 1" in the void and cluster algorithm.
// Initial binary pattern makes blue noise distributed points.
// Phase 1 makes them be progressive, so any points from 0 to N are blue noise.
// Mitchell's best candidate algorithm makes progressive blue noise so can be used instead of those 2 steps.
// https://blog.demofox.org/2017/10/20/generating-blue-noise-sample-points-with-mitchells-best-candidate-algorithm/
static void MitchellsBestCandidate(std::vector<bool>& binaryPattern, std::vector<size_t>& ranks, size_t width, std::mt19937 rng)
{
	ScopedTimer timer("Mitchells Best Candidate", false);

	std::uniform_int_distribution<size_t> dist(0, width * width - 1);

	binaryPattern.resize(width * width, false);
	ranks.resize(width * width, ~size_t(0));

	static const size_t gridCellCount = 32;
	TPointGrid grid(gridCellCount * gridCellCount);
	const size_t gridCellSize = width / gridCellCount;

	size_t ones = size_t(float(width * width) * 0.1f);
	for (size_t i = 0; i < ones; ++i)
	{
		printf("\r%i%%", int(100.0f * float(i) / float(ones - 1)));

		// we scale up the candidates each iteration like in the paper, to keep frequency behavior consistent
		size_t numCandidates = i + 1;

		// keep the candidate that is farthest from the closest existing point
		float bestDistanceSq = 0.0f;
		Point best;
		for (size_t candidate = 0; candidate < numCandidates; ++candidate)
		{
			size_t index = dist(rng);
			Point c;
			c.x = index % width;
			c.y = index / width;

			float minDistSq = DistanceSqToClosestPoint(grid, gridCellCount, gridCellSize, c, width);

			if (minDistSq > bestDistanceSq)
			{
				bestDistanceSq = minDistSq;
				best = c;
			}
		}

		// take the best candidate
		binaryPattern[best.y * width + best.x] = true;
		ranks[best.y * width + best.x] = i;
		AddPointToPointGrid(grid, gridCellCount, gridCellSize, best);
	}
	printf("\n");
}


static size_t VanDerCorput(int base, int index)
{
	double denominator = base, res = 0.0;
	while (index > 0)
	{
		res += (index % base) / denominator;
		index /= base;
		denominator *= base;
	}
	return size_t(res * 64.0f);
}

static size_t Roberts1(uint64_t index)
{
	return size_t(index * 0xC13FA9A902A6328FULL >> 58);
}

static size_t Roberts2(uint64_t index)
{
	return size_t(index * 0x91E10DA5C79E7B1DULL >> 58);
}

static void MakeInitialBinaryPattern(std::vector<bool>& binaryPattern, size_t width, const char* baseFileName, std::mt19937& rng)
{
    ScopedTimer timer("Initial Pattern", false);

    std::uniform_int_distribution<size_t> dist(0, width*width);

    std::vector<float> LUT;
    LUT.resize(width*width, 0.0f);

    binaryPattern.resize(width*width, false);

#if RANDOM_INITIAL()
    size_t ones = size_t(float(width*width) * 0.1f); // start 10% of the pixels as white
    for (size_t index = 0; index < ones; ++index)
    {
        size_t pixel = dist(rng);
        binaryPattern[pixel] = true;
        WriteLUTValue(LUT, width, true, int(pixel % width), int(pixel / width));
    }
#elif VDC_INITIAL()
	size_t x, y;
	for (size_t index = 17; index <= 32; ++index)
	{
		x = VanDerCorput(3, index);
		y = VanDerCorput(5, index);
//		x = Roberts1(uint64_t(index));
//		y = Roberts2(uint64_t(index));
		if (x >= y && x < 63 - y) // bottom edge
		{
			for (size_t xx = 0; xx < 256; xx+=64)
			{
				binaryPattern[xx + x + y * 256] = true;
				binaryPattern[xx + x + (y + 64) * 256] = true;
				WriteLUTValue(LUT, width, true, xx + x, y);
				WriteLUTValue(LUT, width, true, xx + x, y + 64);
			}
		}
		else if (x < y && x >= 63 - y) // top edge
		{
			for (size_t xx = 0; xx < 256; xx += 64)
			{
				binaryPattern[xx + x + y * 256] = true;
				binaryPattern[xx + x + (y + 192) * 256] = true;
				WriteLUTValue(LUT, width, true, xx + x, y);
				WriteLUTValue(LUT, width, true, xx + x, y + 192);
			}
		}
		else if (x < y && x < 63 - y) // left edge
		{
			for (size_t yy = 0; yy < 256; yy += 64)
			{
				binaryPattern[x + (yy + y) * 256] = true;
				binaryPattern[64 + x + (yy + y) * 256] = true;
				WriteLUTValue(LUT, width, true, x, yy + y);
				WriteLUTValue(LUT, width, true, x + 64, yy + y);
			}
		}
		else // right edge
		{
			for (size_t yy = 0; yy < 256; yy += 64)
			{
				binaryPattern[x + (yy + y) * 256] = true;
				binaryPattern[192 + x + (yy + y) * 256] = true;
				WriteLUTValue(LUT, width, true, x, yy + y);
				WriteLUTValue(LUT, width, true, x + 192, yy + y);
			}
		}

		x = VanDerCorput(5, index + 23);
		y = VanDerCorput(7, index + 23);
//		x = VanDerCorput(3, index + 24);
//		y = VanDerCorput(5, index + 24);
//		x = Roberts1(uint64_t(index));
//		y = Roberts2(uint64_t(index));
//		x = Roberts1(uint64_t(index + 32));
//		y = Roberts2(uint64_t(index + 32));
		if (x >= y && x < 63 - y) // bottom edge
		{
			for (size_t xx = 0; xx < 256; xx += 64)
			{
				binaryPattern[xx + x + (y + 128) * 256] = true;
				binaryPattern[xx + x + (y + 192) * 256] = true;
				WriteLUTValue(LUT, width, true, xx + x, y + 128);
				WriteLUTValue(LUT, width, true, xx + x, y + 192);
			}
		}
		else if (x < y && x >= 63 - y) // top edge
		{
			for (size_t xx = 0; xx < 256; xx += 64)
			{
				binaryPattern[xx + x + (y + 64) * 256] = true;
				binaryPattern[xx + x + (y + 128) * 256] = true;
				WriteLUTValue(LUT, width, true, xx + x, y + 64);
				WriteLUTValue(LUT, width, true, xx + x, y + 128);
			}
		}
		else if (x < y && x < 63 - y) // left edge
		{
			for (size_t yy = 0; yy < 256; yy += 64)
			{
				binaryPattern[128 + x + (yy + y) * 256] = true;
				binaryPattern[192 + x + (yy + y) * 256] = true;
				WriteLUTValue(LUT, width, true, x + 128, yy + y);
				WriteLUTValue(LUT, width, true, x + 192, yy + y);
			}
		}
		else // right edge
		{
			for (size_t yy = 0; yy < 256; yy += 64)
			{
				binaryPattern[64 + x + (yy + y) * 256] = true;
				binaryPattern[128 + x + (yy + y) * 256] = true;
				WriteLUTValue(LUT, width, true, x + 64, yy + y);
				WriteLUTValue(LUT, width, true, x + 128, yy + y);
			}
		}
	}
#else
std::vector<size_t> rank1(64 * 64, ~size_t(0));
std::vector<size_t> rank2(64 * 64, ~size_t(0));
std::vector<bool> ibp1, ibp2;
MitchellsBestCandidate(ibp1, rank1, 64, rng);
for (size_t x = 0; x < 64; x++)
{
	for (size_t y = 0; y < 64; y++)
	{
		if (rank1[x + y * 64] > 23)
			continue;
		if (x >= y && x < 63 - y) // bottom edge
		{
			for (size_t xx = 0; xx < 256; xx += 64)
			{
				binaryPattern[xx + x + y * 256] = true;
				binaryPattern[xx + x + (y + 64) * 256] = true;
				WriteLUTValue(LUT, width, true, xx + x, y);
				WriteLUTValue(LUT, width, true, xx + x, y + 64);
			}
		}
		else if (x < y && x >= 63 - y) // top edge
		{
			for (size_t xx = 0; xx < 256; xx += 64)
			{
				binaryPattern[xx + x + y * 256] = true;
				binaryPattern[xx + x + (y + 192) * 256] = true;
				WriteLUTValue(LUT, width, true, xx + x, y);
				WriteLUTValue(LUT, width, true, xx + x, y + 192);
			}
		}
		else if (x < y && x < 63 - y) // left edge
		{
			for (size_t yy = 0; yy < 256; yy += 64)
			{
				binaryPattern[x + (yy + y) * 256] = true;
				binaryPattern[64 + x + (yy + y) * 256] = true;
				WriteLUTValue(LUT, width, true, x, yy + y);
				WriteLUTValue(LUT, width, true, x + 64, yy + y);
			}
		}
		else // right edge
		{
			for (size_t yy = 0; yy < 256; yy += 64)
			{
				binaryPattern[x + (yy + y) * 256] = true;
				binaryPattern[192 + x + (yy + y) * 256] = true;
				WriteLUTValue(LUT, width, true, x, yy + y);
				WriteLUTValue(LUT, width, true, x + 192, yy + y);
			}
		}
	}
}
MitchellsBestCandidate(ibp2, rank2, 64, rng);
for (size_t x = 0; x < 64; x++)
{
	for (size_t y = 0; y < 64; y++)
	{
		if (rank2[x + y * 64] > 42)
			continue;
		if (x >= y && x < 63 - y) // bottom edge
		{
			for (size_t xx = 0; xx < 256; xx += 64)
			{
				binaryPattern[xx + x + (y + 128) * 256] = true;
				binaryPattern[xx + x + (y + 192) * 256] = true;
				WriteLUTValue(LUT, width, true, xx + x, y + 128);
				WriteLUTValue(LUT, width, true, xx + x, y + 192);
			}
		}
		else if (x < y && x >= 63 - y) // top edge
		{
			for (size_t xx = 0; xx < 256; xx += 64)
			{
				binaryPattern[xx + x + (y + 64) * 256] = true;
				binaryPattern[xx + x + (y + 128) * 256] = true;
				WriteLUTValue(LUT, width, true, xx + x, y + 64);
				WriteLUTValue(LUT, width, true, xx + x, y + 128);
			}
		}
		else if (x < y && x < 63 - y) // left edge
		{
			for (size_t yy = 0; yy < 256; yy += 64)
			{
				binaryPattern[128 + x + (yy + y) * 256] = true;
				binaryPattern[192 + x + (yy + y) * 256] = true;
				WriteLUTValue(LUT, width, true, x + 128, yy + y);
				WriteLUTValue(LUT, width, true, x + 192, yy + y);
			}
		}
		else // right edge
		{
			for (size_t yy = 0; yy < 256; yy += 64)
			{
				binaryPattern[64 + x + (yy + y) * 256] = true;
				binaryPattern[128 + x + (yy + y) * 256] = true;
				WriteLUTValue(LUT, width, true, x + 64, yy + y);
				WriteLUTValue(LUT, width, true, x + 128, yy + y);
			}
		}
	}
}

#endif
    int iterationCount = 0;
    while (1)
    {
        printf("\r%i iterations", iterationCount);
        iterationCount++;

        // find the location of the tightest cluster
        int tightestClusterX = -1;
        int tightestClusterY = -1;
        FindTightestClusterLUT(LUT, binaryPattern, width, tightestClusterX, tightestClusterY, rng);

        // remove the 1 from the tightest cluster
        binaryPattern[tightestClusterY*width + tightestClusterX] = false;
        WriteLUTValue(LUT, width, false, tightestClusterX, tightestClusterY);

        // find the largest void
        int largestVoidX = -1;
        int largestVoidY = -1;
        FindLargestVoidLUT(LUT, binaryPattern, width, largestVoidX, largestVoidY, rng);

        // put the 1 in the largest void
        binaryPattern[largestVoidY*width + largestVoidX] = true;
        WriteLUTValue(LUT, width, true, largestVoidX, largestVoidY);

        #if SAVE_VOIDCLUSTER_INITIALBP()
        // save the binary pattern out for debug purposes
        SaveBinaryPattern(binaryPattern, width, baseFileName, iterationCount, tightestClusterX, tightestClusterY, largestVoidX, largestVoidY);
        #endif

        // exit condition. the pattern is stable
        if (tightestClusterX == largestVoidX && tightestClusterY == largestVoidY)
            break;
    }
    printf("\n");
}

// Phase 1: Start with initial binary pattern and remove the tightest cluster until there are none left, entering ranks for those pixels
static void Phase1(std::vector<bool>& binaryPattern, std::vector<float>& LUT, std::vector<size_t>& ranks, size_t width, std::mt19937& rng, const char* baseFileName)
{
    ScopedTimer timer("Phase 1", false);

    // count how many ones there are
    size_t ones = 0;
    for (bool b : binaryPattern)
    {
        if (b)
            ones++;
    }
    size_t startingOnes = ones;

    // remove the tightest cluster repeatedly
    while (ones > 0)
    {
        printf("\r%i%%", int(100.0f * (1.0f - float(ones) / float(startingOnes))));

        int bestX, bestY;
        FindTightestClusterLUT(LUT, binaryPattern, width, bestX, bestY, rng);
        binaryPattern[bestY * width + bestX] = false;
        WriteLUTValue(LUT, width, false, bestX, bestY);
        ones--;
        ranks[bestY*width + bestX] = ones;

        #if SAVE_VOIDCLUSTER_PHASE1()
        // save the binary pattern out for debug purposes
        SaveBinaryPattern(binaryPattern, width, baseFileName, int(startingOnes - ones), bestX, bestY, -1, -1);
        #endif
    }
    printf("\n");
}


// Phase 2: Start with initial binary pattern and add points to the largest void until half the pixels are white, entering ranks for those pixels
static void Phase2(std::vector<bool>& binaryPattern, std::vector<float>& LUT, std::vector<size_t>& ranks, size_t width, std::mt19937& rng)
{
    ScopedTimer timer("Phase 2", false);

    // count how many ones there are
    size_t ones = 0;
    for (bool b : binaryPattern)
    {
        if (b)
            ones++;
    }
    size_t startingOnes = ones;
    size_t onesToDo = (width*width / 2) - startingOnes;

    // add to the largest void repeatedly
    while (ones <= (width*width/2))
    {
        size_t onesDone = ones - startingOnes;
        printf("\r%i%%", int(100.0f * float(onesDone) / float(onesToDo)));

        int bestX, bestY;
        FindLargestVoidLUT(LUT, binaryPattern, width, bestX, bestY, rng);
        binaryPattern[bestY * width + bestX] = true;
        WriteLUTValue(LUT, width, true, bestX, bestY);
        ranks[bestY*width + bestX] = ones;
        ones++;
    }
    printf("\n");
}

// Phase 3: Continue with the last binary pattern, repeatedly find the tightest cluster of 0s and insert a 1 into them
static void Phase3(std::vector<bool>& binaryPattern, std::vector<float>& LUT, std::vector<size_t>& ranks, size_t width, std::mt19937& rng)
{
    ScopedTimer timer("Phase 3", false);

    // count how many ones there are
    size_t ones = 0;
    for (bool b : binaryPattern)
    {
        if (b)
            ones++;
    }
    size_t startingOnes = ones;
    size_t onesToDo = (width*width) - startingOnes;

    // add 1 to the largest cluster of 0's repeatedly
    int bestX, bestY;
    while (FindLargestVoidLUT(LUT, binaryPattern, width, bestX, bestY, rng))
    {
        size_t onesDone = ones - startingOnes;
        printf("\r%i%%", int(100.0f * float(onesDone) / float(onesToDo)));

        WriteLUTValue(LUT, width, true, bestX, bestY);
        binaryPattern[bestY * width + bestX] = true;
        ranks[bestY*width + bestX] = ones;
        ones++;
    }
    printf("\n");
}

void TestMask(const std::vector<uint8_t>& noise, size_t noiseSize, const char* baseFileName)
{
	std::vector<uint8_t> thresholdImage(noise.size());

	for (size_t testIndex = 0; testIndex < THRESHOLD_SAMPLES(); ++testIndex)
	{
		float percent = float(testIndex) / float(THRESHOLD_SAMPLES() - 1);
		uint8_t thresholdValue = FromFloat<uint8_t>(percent);
		if (thresholdValue == 0)
			thresholdValue = 1;
		else if (thresholdValue == 255)
			thresholdValue = 254;

		for (size_t pixelIndex = 0, pixelCount = noise.size(); pixelIndex < pixelCount; ++pixelIndex)
			thresholdImage[pixelIndex] = noise[pixelIndex] > thresholdValue ? 255 : 0;

		std::vector<uint8_t> thresholdImageDFT;
		DFT(thresholdImage, thresholdImageDFT, noiseSize);

		std::vector<uint8_t> noiseAndDFT;
		size_t noiseAndDFT_width = 0;
		size_t noiseAndDFT_height = 0;
		AppendImageHorizontal(thresholdImage, noiseSize, noiseSize, thresholdImageDFT, noiseSize, noiseSize, noiseAndDFT, noiseAndDFT_width, noiseAndDFT_height);

		char fileName[256];
		sprintf(fileName, "%s_%u.png", baseFileName, thresholdValue);
		stbi_write_png(fileName, int(noiseAndDFT_width), int(noiseAndDFT_height), 1, noiseAndDFT.data(), 0);
	}
}

void TestNoise(const std::vector<uint8_t>& noise, size_t noiseSize, const char* baseFileName)
{
	char fileName[256];
	sprintf(fileName, "%s.histogram.csv", baseFileName);

	WriteHistogram(noise, fileName);
	std::vector<uint8_t> noiseDFT;
	DFT(noise, noiseDFT, noiseSize);

	//std::vector<uint8_t> noiseAndDFT;
	//size_t noiseAndDFT_width = 0;
	//size_t noiseAndDFT_height = 0;
	//AppendImageHorizontal(noise, noiseSize, noiseSize, noiseDFT, noiseSize, noiseSize, noiseAndDFT, noiseAndDFT_width, noiseAndDFT_height);

	sprintf(fileName, "%s.png", baseFileName);
	stbi_write_png(fileName, int(noiseSize), int(noiseSize), 1, noise.data(), 0);

	sprintf(fileName, "%s_hist.png", baseFileName);
	stbi_write_png(fileName, int(noiseSize), int(noiseSize), 1, noiseDFT.data(), 0);

#if TEST_MASK()
	TestMask(noise, noiseSize, baseFileName);
#endif
}


void GenerateBN_Void_Cluster(std::vector<uint8_t>& blueNoise, size_t width, bool useMitchellsBestCandidate, const char* baseFileName)
{
	std::mt19937 rng(GetRNGSeed());
	std::vector<size_t> ranks(width * width, ~size_t(0));

	std::vector<bool> initialBinaryPattern;
	std::vector<bool> binaryPattern;
	std::vector<float> initialLUT;
	std::vector<float> LUT;

	if (!useMitchellsBestCandidate)
	{
		// make the initial binary pattern and initial LUT
		MakeInitialBinaryPattern(initialBinaryPattern, width, baseFileName, rng);
		MakeLUT(initialBinaryPattern, initialLUT, width, true);

		// Phase 1: Start with initial binary pattern and remove the tightest cluster until there are none left, entering ranks for those pixels
		binaryPattern = initialBinaryPattern;
		LUT = initialLUT;
		Phase1(binaryPattern, LUT, ranks, width, rng, baseFileName);
	}
	else
	{
		// replace initial binary pattern and phase 1 with Mitchell's best candidate algorithm, and then making the LUT
		MitchellsBestCandidate(initialBinaryPattern, ranks, width, rng);
		MakeLUT(initialBinaryPattern, initialLUT, width, true);

		//SaveBinaryPattern(initialBinaryPattern, width, "out/_blah", 0, -1, -1, -1, -1);
	}

	// Phase 2: Start with initial binary pattern and add points to the largest void until half the pixels are white, entering ranks for those pixels
	binaryPattern = initialBinaryPattern;
	LUT = initialLUT;
	Phase2(binaryPattern, LUT, ranks, width, rng);

	// Phase 3: Continue with the last binary pattern, repeatedly find the tightest cluster of 0s and insert a 1 into them
	// Note: we do need to re-make the LUT, because we are writing 0s instead of 1s
	MakeLUT(binaryPattern, LUT, width, false);
	Phase3(binaryPattern, LUT, ranks, width, rng);

	// convert to U8
	{
		ScopedTimer timer("Converting to U8", false);
		blueNoise.resize(width * width);
		for (size_t index = 0; index < width * width; ++index)
			blueNoise[index] = uint8_t(ranks[index] * 256 / (width * width));
	}
	TestNoise(blueNoise, width, baseFileName);
}
