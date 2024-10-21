#include <random>
#include<iostream>
#include <spdlog/spdlog.h>
#include<vector>

#include "Labs/1-Drawing2D/tasks.h"

using VCX::Labs::Common::ImageRGB;
using namespace std;



namespace VCX::Labs::Drawing2D {
    /******************* 1.Image Dithering *****************/
    void DitheringThreshold(
        ImageRGB &       output,
        ImageRGB const & input) {
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color = input.At(x, y);
                output.At(x, y) = {
                    color.r > 0.5 ? 1 : 0,
                    color.g > 0.5 ? 1 : 0,
                    color.b > 0.5 ? 1 : 0,
                };
            }
    }

    void DitheringRandomUniform(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        for (size_t x = 0; x < input.GetSizeX(); ++x)
            for (size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color = input.At(x, y);
                float     ran   = float(rand() % 10000) / 10000 - 0.5;
                // cout << ran << endl;
                color.r += ran;
                color.g += ran;
                color.b += ran;
                output.At(x, y) = {
                    color.r > 0.5 ? 1 : 0,
                    color.g > 0.5 ? 1 : 0,
                    color.b > 0.5 ? 1 : 0,
                };
            }
    }

    void DitheringRandomBlueNoise(
        ImageRGB &       output,
        ImageRGB const & input,
        ImageRGB const & noise) {
        // your code here:
        for (size_t x = 0; x < input.GetSizeX(); ++x)
            for (size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color = input.At(x, y);
                glm::vec3 n     = noise.At(x, y);
                color.r += n.r;
                color.g += n.g;
                color.b += n.b;
                output.At(x, y) = {
                    color.r - 0.5 > 0.5 ? 1 : 0,
                    color.g - 0.5 > 0.5 ? 1 : 0,
                    color.b - 0.5 > 0.5 ? 1 : 0,
                };
            }
    }

    void DitheringOrdered(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        int lst[9][2] = {
            { 1, 1 },
            { 0, 1 },
            { 1, 2 },
            { 2, 1 },
            { 2, 0 },
            { 0, 2 },
            { 0, 0 },
            { 2, 2 },
            { 1, 0 }
        };
        // your code here:
        for (size_t x = 0; x < input.GetSizeX(); ++x)
            for (size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color      = input.At(x, y);
                int       degree     = int(color.r * 9);
                int       new3[3][3] = { 0 };
                for (int i = 0; i <= degree; i++) {
                    int x_       = lst[i][0];
                    int y_       = lst[i][1];
                    new3[x_][y_] = 1;
                }

                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++) {
                        output.At(3 * x + i, 3 * y + j) = {
                            new3[i][j],
                            new3[i][j],
                            new3[i][j]
                        };
                    }
            }
    }

    void DitheringErrorDiffuse(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        glm::vec3 arr[180][215] ;
        for (int x = 0; x < input.GetSizeX(); ++x)
            for (int y = 0; y < input.GetSizeY(); ++y) {
                // output.At(x, y) = input.At(x, y);
                arr[x][y] = input.At(x, y);
            }
        for (std::size_t y = 0; y < input.GetSizeY(); ++y)
            for (std::size_t x = 0; x < input.GetSizeX(); ++x) {
                glm::vec3 color = arr[x][y];
                arr[x][y] = {
                    color.r > 0.5 ? 1.0 : 0.0,
                    color.g > 0.5 ? 1.0 : 0.0,
                    color.b > 0.5 ? 1.0 : 0.0,
                };
                glm::vec3 color_ = arr[x][y];
                float     error  = color.r - color_.r;

                if (x < input.GetSizeX() - 1) {
                    glm::vec3 color_r = arr[x+1][y];
                    color_r += error * 7.0f / 16.0f;
                    arr[x + 1][y] = color_r;
                    ;

                    if (y < input.GetSizeY() - 1) {
                        glm::vec3 color_rd = arr[x+1][y+1];
                        color_rd += error / 16.0f;
                        arr[x + 1][y + 1] = color_rd;

                        // glm::vec3 color_rd = output.At(x + 1, y + 1);
                    }
                }
                if (y < input.GetSizeY() - 1) {
                    glm::vec3 color_d = arr[x][y+1];
                    color_d += error * 5.0f / 16.0f;
                    arr[x][y+1]=color_d;

                    if (x > 0) {
                        glm::vec3 color_ld = arr[x-1][y+1];
                        color_ld += error * 3.0f / 16.0f;
                        arr[x-1][y+1]=color_ld;
                    }
                }
            }
        for (int x = 0; x < input.GetSizeX();x++){
            for (int y = 0; y < input.GetSizeY();y++){
                glm::vec3 color = arr[x][y];
                output.At(x, y) = {
                    color.r > 0.5 ? 1 : 0,
                    color.g > 0.5 ? 1 : 0,
                    color.b > 0.5 ? 1 : 0,
                };
            }
        }
    }

    /******************* 2.Image Filtering *****************/
    void Blur(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        float kernel[3][3] = {
            { 1.0f / 9.0f, 1.0f / 9.0f, 1.0f / 9.0f },
            { 1.0f / 9.0f, 1.0f / 9.0f, 1.0f / 9.0f },
            { 1.0f / 9.0f, 1.0f / 9.0f, 1.0f / 9.0f }
        };
        for (size_t x = 0; x < input.GetSizeX() - 2; x++) {
            for (size_t y = 0; y < input.GetSizeY() - 2; y++) {
                glm::vec3 color[3][3] = {
                    {     input.At(x, y),     input.At(x, y + 1),     input.At(x, y + 2) },
                    { input.At(x + 1, y), input.At(x + 1, y + 1), input.At(x + 1, y + 2) },
                    { input.At(x + 2, y), input.At(x + 2, y + 1), input.At(x + 2, y + 2) }
                };

                glm::vec3 o_color = { 0.0, 0.0, 0.0 };
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        o_color += kernel[i][j] * color[i][j];
                    }
                }
                output.At(x, y) = o_color;
            }
        }
    }

    void Edge(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        float Gx[3][3] = {
            { -1.0, 0.0, 1.0 },
            { -2.0, 0.0, 2.0 },
            { -1.0, 0.0, 1.0 }
        };
        float Gy[3][3] = {
            {  1.0,  2.0,  1.0 },
            {  0.0,  0.0,  0.0 },
            { -1.0, -2.0, -1.0 }
        };
        for (size_t x = 0; x < input.GetSizeX() - 2; x++) {
            for (size_t y = 0; y < input.GetSizeY() - 2; y++) {
                glm::vec3 color[3][3] = {
                    {     input.At(x, y),     input.At(x, y + 1),     input.At(x, y + 2) },
                    { input.At(x + 1, y), input.At(x + 1, y + 1), input.At(x + 1, y + 2) },
                    { input.At(x + 2, y), input.At(x + 2, y + 1), input.At(x + 2, y + 2) }
                };

                glm::vec3 o_colorx = { 0.0, 0.0, 0.0 };
                glm::vec3 o_colory = { 0.0, 0.0, 0.0 };
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        o_colorx += Gx[i][j] * color[i][j];
                        o_colory += Gy[i][j] * color[i][j];
                    }
                }
                glm::vec3 o_color = sqrt(o_colory * o_colory + o_colorx * o_colorx);
                output.At(x, y)   = o_color;
            }
        }
    }

    /******************* 3. Image Inpainting *****************/
    void Inpainting(
        ImageRGB &         output,
        ImageRGB const &   inputBack,
        ImageRGB const &   inputFront,
        const glm::ivec2 & offset) {
        output             = inputBack;
        std::size_t width  = inputFront.GetSizeX();
        std::size_t height = inputFront.GetSizeY();
        glm::vec3 * g      = new glm::vec3[width * height];
        memset(g, 0, sizeof(glm::vec3) * width * height);
        // set boundary condition
        for (std::size_t y = 0; y < height; ++y) {
            // set boundary for (0, y), your code: g[y * width] = ?
            g[y * width] = inputBack.At(offset.x, offset.y + y) - inputFront.At(0, y);
            // set boundary for (width - 1, y), your code: g[y * width + width - 1] = ?
            g[y * width + width - 1] = inputBack.At(offset.x + width - 1, offset.y + y) - inputFront.At(width - 1, y);
        }
        for (std::size_t x = 0; x < width; ++x) {
            // set boundary for (x, 0), your code: g[x] = ?
            g[x] = inputBack.At(offset.x + x, offset.y) - inputFront.At(x, 0);
            // set boundary for (x, height - 1), your code: g[(height - 1) * width + x] = ?
            g[(height - 1) * width + x] = inputBack.At(offset.x + x, offset.y + height - 1) - inputFront.At(x, height - 1);
        }

        // Jacobi iteration, solve Ag = b
        for (int iter = 0; iter < 8000; ++iter) {
            for (std::size_t y = 1; y < height - 1; ++y)
                for (std::size_t x = 1; x < width - 1; ++x) {
                    g[y * width + x] = (g[(y - 1) * width + x] + g[(y + 1) * width + x] + g[y * width + x - 1] + g[y * width + x + 1]);
                    g[y * width + x] = g[y * width + x] * glm::vec3(0.25);
                }
        }

        for (std::size_t y = 0; y < inputFront.GetSizeY(); ++y)
            for (std::size_t x = 0; x < inputFront.GetSizeX(); ++x) {
                glm::vec3 color = g[y * width + x] + inputFront.At(x, y);
                output.At(x + offset.x, y + offset.y) = color;
            }
        delete[] g;
    }

    /******************* 4. Line Drawing *****************/
    void DrawLine(
        ImageRGB &       canvas,
        glm::vec3 const  color,
        glm::ivec2 const p0,
        glm::ivec2 const p1) {
        // your code here:
        int x0 = p0[0], y0 = p0[1];
        int x1 = p1[0], y1 = p1[1];
        int dx = x1 - x0;
        int dy = y1 - y0;
        int ux = dx > 0 ? 1 : -1;
        int uy = dy > 0 ? 1 : -1;
        int x = x0, y = y0, eps = 0;
        dx = abs(dx);
        dy = abs(dy);
        if (dx > dy) {
            for (x = x0; x != x1; x += ux) {
                canvas.At(x, y) = color;
                eps += dy;
                if (2 * eps >= dx) {
                    y += uy;
                    eps -= dx;
                }
            }
        } else {
            for (y = y0; y != y1; y += uy) {
                canvas.At(x, y) = color;
                eps += dx;
                if (2 * eps >= dy) {
                    x += ux;
                    eps -= dy;
                }
            }
        }
    }

    /******************* 5. Triangle Drawing *****************/
    void DrawTriangleFilled(
        ImageRGB &       canvas,
        glm::vec3 const  color,
        glm::ivec2 const p0,
        glm::ivec2 const p1,
        glm::ivec2 const p2) {
        // your code here:
        float                                   EPS    = 0.0001;
        std::vector<std::pair<glm::ivec2, int>> points = {
            { p0, p0.y },
            { p1, p1.y },
            { p2, p2.y }
        };

        // 按照 y 值进行排序
        std::sort(points.begin(), points.end(), [](const auto & a, const auto & b) {
            return a.second < b.second; // 升序排序
        });

        // points 现在是按 y 值从小到大排序的
        glm::ivec2 first  = points[0].first; // y 值最小的点1
        glm::ivec2 second = points[1].first; // y 值第二大的点2
        glm::ivec2 third  = points[2].first; // y 值最大的点3
        
        int        y;

        float d12_x = second[0] - first[0];
        float d12_y = second[1] - first[1];
        float d23_x = third[0] - second[0];
        float d23_y = third[1] - second[1];
        float d13_x = third[0] - first[0];
        float d13_y = third[1] - first[1];
        float d12   = d12_x / (d12_y + EPS);//三条边的斜率
        float d23   = d23_x / (d23_y + EPS);
        float d13   = d13_x / (d13_y + EPS);
        float xl = first[0], xr = first[0];
        for (y = first[1]; y < second[1]; y++) {//y值从first[1]到second[1]
            for (int x = xl; x <= xr; x++) {
                canvas.At(x, y) = color;
            }
            if (d12 > d13) {
                xl += d13;
                xr += d12;
            } else {
                xl += d12;
                xr += d13;
            }
        }
        xl = third[0];
        xr = third[0];
        for (y = third[1]; y >= second[1]; y--) {//y值从third[1]到second[1]
            for (int x = xl; x <= xr; x++) {
                canvas.At(x, y) = color;
            }
            if (d13 > d23) {
                xl -= d13;
                xr -= d23;
            } else {
                xl -= d23;
                xr -= d13;
            }
        }
    
    }

    /******************* 6. Image Supersampling *****************/
    void Supersample(
        ImageRGB &       output,
        ImageRGB const & input,
        int              rate) {
        // your code here:
        int in_width   = input.GetSizeX();
        int in_height  = input.GetSizeY();
        int out_width  = in_width / rate;
        int out_height = in_height / rate;

        for (size_t x = 0; x < out_width; x++) {
            for (size_t y = 0; y < out_height; y++) {
                glm::vec3 color = { 0.0, 0.0, 0.0 };
                for (int i = 0; i < rate; i++) {
                    for (int j = 0; j < rate; j++) {
                        color += input.At(x * rate + i, y * rate + j);
                    }
                }
                color = color / float(rate * rate);
                // input.At(x * rate, y * rate) = color;
            }
        }
        float dx = float(out_width) / float(320);
        float dy = float(out_height) / float(320);
        for (float x = 0.0; x < 320; x++) {
            for (float y = 0.0; y < 320; y++) {
                glm::vec3 color = { 0.0, 0.0, 0.0 };
                for (int i = 0; i < rate; i++) {
                    for (int j = 0; j < rate; j++) {
                        color += input.At(int(x * dx) * rate+i, int(y * dy) * rate+j);
                    }
                }
                color           = color / float(rate * rate);
                output.At(x, y) = color;
            }
        }
    }

    /******************* 7. Bezier Curve *****************/
    // Note: Please finish the function [DrawLine] before trying this part.
    glm::vec2 CalculateBezierPoint(
        std::span<glm::vec2> points,
        float const          t) {
        // your code here:
        int len = points.size();
        vector<glm::vec2> points_(len);
        // copy(points.begin(), points.end(), points_.begin());
        for (int i = 0; i < len; i++) {
            points_[i] = points[i];
        }
        for (int n = len; n >= 2; n--) {
            for (int i = 0; i < n - 1; i++) {
                glm::vec2 temp = points_[i] * glm::vec2(1.0 - t) + points_[i + 1] * glm::vec2(t);
                points_[i]     = temp;
            }
        }
        return points_[0];
    }
} // namespace VCX::Labs::Drawing2D