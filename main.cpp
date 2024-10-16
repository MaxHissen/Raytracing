// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
// Utilities for the Assignment
#include "raster.h"
#include <gif.h>
#include <fstream>
#include <Eigen/Geometry>
// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"
using namespace std;
using namespace Eigen;

//Image height
const int H = 480;

//Camera settings
const double near_plane = 1.5; //AKA focal length
const double far_plane = near_plane * 100;
const double field_of_view = 0.7854; //45 degrees
const double aspect_ratio = 1.5;
const bool is_perspective = false;
const Vector3d camera_position(0, 0, 3);
const Vector3d camera_gaze(0, 0, -1);
const Vector3d camera_top(0, 1, 0);

//Object
const string data_dir = DATA_DIR;
const string mesh_filename(data_dir + "bunny.off");
MatrixXd vertices;
MatrixXi facets;

//Material for the object
const Vector3d obj_diffuse_color(0.5, 0.5, 0.5);
const Vector3d obj_specular_color(0.2, 0.2, 0.2);
const double obj_specular_exponent = 256.0;

//Lights
std::vector<Vector3d> light_positions;
std::vector<Vector3d> light_colors;
//Ambient light
const Vector3d ambient_light(0.3, 0.3, 0.3);

//Fills the different arrays
void setup_scene(){

    std::ifstream in(mesh_filename);
    if (!in.good()){
        std::cerr << "Invalid file " << mesh_filename << std::endl;
        exit(1);
    }
    std::string token;
    in >> token;
    int nv, nf, ne;
    in >> nv >> nf >> ne;
    vertices.resize(nv, 3);
    facets.resize(nf, 3);
    for (int i = 0; i < nv; ++i){
        in >> vertices(i, 0) >> vertices(i, 1) >> vertices(i, 2);
    }
    for (int i = 0; i < nf; ++i){
        int s;
        in >> s >> facets(i, 0) >> facets(i, 1) >> facets(i, 2);
        assert(s == 3);
    }

    //Lights
    light_positions.emplace_back(8, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(6, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(4, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(2, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(0, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(-2, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(-4, 8, 0);
    light_colors.emplace_back(16, 16, 16);
}

void build_uniform(UniformAttributes &uniform){
    // setup uniform

    // setup camera, compute w, u, v
    Vector3d w = -camera_gaze / camera_gaze.norm();
    Vector3d u = camera_top.cross(w) / (camera_top.cross(w)).norm();
    Vector3d v = w.cross(u);

    Vector4d w4 ( w(0), w(1), w(2), 0 );
    Vector4d u4 ( u(0), u(1), u(2), 0);
    Vector4d v4 ( v(0), v(1), v(2), 0);
    Vector4d e4 ( camera_position(0), camera_position(1), camera_position(2), 1);

    // compute the camera transformation
    MatrixXd Mtrx_Camera;
    Mtrx_Camera.resize(4, 4);
    Mtrx_Camera.col(0) = u4;
    Mtrx_Camera.col(1) = v4;
    Mtrx_Camera.col(2) = w4;
    Mtrx_Camera.col(3) = e4;
    Mtrx_Camera = Mtrx_Camera.inverse();

    // setup projection matrix
    double t = tan(field_of_view / 2.0f) * near_plane;
    double r = t * aspect_ratio; 
    double l = -r;
    double b = -t;
    double n = -near_plane;
    double f = -far_plane;

    MatrixXd Orth_Matrix;
    Orth_Matrix.resize(4,4);

    Orth_Matrix <<  Vector4d(2/(r-l), 0, 0, -(r+l)/(r-l)), 
                    Vector4d(0, 2/(t-b), 0, -(t+b) / (t-b)), 
                    Vector4d(0, 0, 2 / (n-f), -(n+f)/(n-f)), 
                    Vector4d(0, 0, 0, 1);
    Orth_Matrix.transposeInPlace();

    

    Matrix4d P;
    P <<    Vector4d(n, 0, 0, 0), 
            Vector4d(0, n, 0, 0),
            Vector4d(0, 0, n+f, -f*n),
            Vector4d(0, 0, 1, 0);
    P.transposeInPlace();

    if (is_perspective){
        // setup prespective camera
        uniform.view = (Orth_Matrix * P * Mtrx_Camera);
    }else{
        uniform.view = (Orth_Matrix * Mtrx_Camera);
    }
}

void simple_render(Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer){
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;

    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        VertexAttributes transformed_VA;
        transformed_VA.position = uniform.view * va.position; 
        return transformed_VA;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        return FragmentAttributes(1, 0, 0);
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        return FrameBufferAttributes(fa.color[0] * 255, fa.color[1] * 255, fa.color[2] * 255, fa.color[3] * 255);
    };

    std::vector<VertexAttributes> vertex_attributes;

        for(int i = 0; i < facets.rows(); i++) {
        Vector3i triangle = facets.row(i);

        for(int j = 0; j < 3; j++) {

            Vector3d v = vertices.row(triangle[j]);

            vertex_attributes.push_back(VertexAttributes(v[0], v[1], v[2]));
        }
    }
    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

Matrix4d compute_rotation(const double alpha){
    // Compute the rotation matrix of angle alpha on the y axis around the object barycenter
    Matrix4d res;
    res <<  Vector4d(cos(alpha), 0, -sin(alpha), 0), 
            Vector4d(0, 1, 0, 0),
            Vector4d(sin(alpha), 0, cos(alpha), 0), 
            Vector4d(0, 0, 0, 1);
    
    return res;
}

void wireframe_render(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer){
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;

    Matrix4d trafo = compute_rotation(alpha);
    uniform.transform = trafo;

    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        // fill the shader
        VertexAttributes transformed;
        transformed.position = uniform.view * (uniform.transform * va.position);
        return transformed;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        // fill the shader
        return FragmentAttributes(1, 0, 0);
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        // fill the shader
        return FrameBufferAttributes(fa.color[0] * 255, fa.color[1] * 255, fa.color[2] * 255, fa.color[3] * 255);
    };

    std::vector<VertexAttributes> vertex_attributes;

    // generate the vertex attributes for the edges and rasterize the lines
    // use the transformation matrix

    for (int i = 0; i < facets.rows(); i++) {
        Vector3i r = facets.row(i);
        for (int j = 0; j < 3; j ++) {
            VertexAttributes v(vertices(r(j), 0), vertices(r(j), 1), vertices(r(j), 2));
            vertex_attributes.push_back(v);
            int k = (j == 2) ? 0 : j + 1;
            VertexAttributes b(vertices(r(k), 0), vertices(r(k), 1), vertices(r(k), 2));
            vertex_attributes.push_back(b);
        }
    }

    rasterize_lines(program, uniform, vertex_attributes, 0.5, frameBuffer);
}

void get_shading_program(Program &program){
   program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        // transform the position and the normal
        // compute the correct lighting
        VertexAttributes new_va;

        new_va.position = uniform.view * (uniform.transform * va.position);
        new_va.normal = uniform.view * va.normal;

        const Vector3d pos (new_va.position(0), new_va.position(1), new_va.position(2));

        Vector4d n_transform = uniform.transform * va.normal;

        const Vector3d norm (n_transform(0), n_transform(1), n_transform(2));

        Vector3d v = (camera_position - pos).normalized();

        Vector3d lights_color(0,0,0);

        for (int i = 0; i < light_positions.size(); i ++) {
            const Vector3d &cur_light_position = light_positions[i];
            const Vector3d &cur_light_color = light_colors[i];

            const Vector3d Li = (cur_light_position - pos).normalized();
            const Vector3d diffuse = obj_diffuse_color * std::max(Li.dot(norm), 0.0);

            Vector3d l = (cur_light_position - pos).normalized();
            Vector3d h = (v + l) / (v + l).norm();
            const Vector3d specular = obj_specular_color * std::pow(h.transpose() * norm, obj_specular_exponent);

            Vector3d D = cur_light_position - pos;
            lights_color += (diffuse + specular).cwiseProduct(cur_light_color) / D.squaredNorm();
        }

        new_va.color = ambient_light + lights_color;
        return new_va;
    }; 

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        // create the correct fragment
 
        FragmentAttributes output(va.color(0), va.color(1), va.color(2));
        if(is_perspective) {
            output.position = va.position;
        }else{
            output.position = va.position.cwiseProduct(Vector4d(1, 1, -1, 1));
        }
        return output;
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        // implement the depth check
        if (fa.position[2] < previous.depth){
            FrameBufferAttributes new_fa(fa.color[0] * 255, fa.color[1] * 255, fa.color[2] * 255, fa.color[3] * 255);
            new_fa.depth = fa.position[2];
            return new_fa;
            
        }else{
            return previous;
        }
    };
}

void flat_shading(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer){
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;
    get_shading_program(program);

    Eigen::Matrix4d trafo = compute_rotation(alpha);
    uniform.transform = trafo;
    std::vector<VertexAttributes> vertex_attributes;
    // compute the normals
    // set material colors
    for (int i = 0; i < facets.rows(); i++) {
        Vector3i r = facets.row(i);

        Vector3d a (vertices(r(0), 0), vertices(r(0), 1), vertices(r(0), 2));
        Vector3d b (vertices(r(1), 0), vertices(r(1), 1), vertices(r(1), 2));
        Vector3d c (vertices(r(2), 0), vertices(r(2), 1), vertices(r(2), 2));

        Vector3d n = ((b - a).cross(c - a)).normalized();
        Vector4d normal (n(0), n(1), n(2), 0);

        for (int j = 0; j < 3; j ++) {
            VertexAttributes v(vertices(r(j), 0), vertices(r(j), 1), vertices(r(j), 2));
            v.normal = normal;
            vertex_attributes.push_back(v);
        }
    }

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

void pv_shading(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;
    get_shading_program(program);

    Eigen::Matrix4d trafo = compute_rotation(alpha);
    uniform.transform = trafo;

    // compute the vertex normals as vertex normal average

    std::vector<VertexAttributes> vertex_attributes;
    Eigen::MatrixXd Vector_Normal;
    Vector_Normal.resize(vertices.rows(), 3);
    Vector_Normal.setZero();
    
    for(int i = 0; i < facets.rows(); i ++) {
        Vector3i r = facets.row(i);
        Vector3d a (vertices(r(0), 0), vertices(r(0), 1), vertices(r(0), 2));
        Vector3d b (vertices(r(1), 0), vertices(r(1), 1), vertices(r(1), 2));
        Vector3d c (vertices(r(2), 0), vertices(r(2), 1), vertices(r(2), 2));

        Vector3d n = ((b - a).cross(c - a)).normalized();

        Vector_Normal.row(r(0)) += n;
        Vector_Normal.row(r(1)) += n;
        Vector_Normal.row(r(1)) += n;
    }

    for (int i = 0; i < facets.rows(); i++) {
        Vector3i r = facets.row(i);
        for (int j = 0; j < 3; j ++) {
            VertexAttributes v(vertices(r(j), 0), vertices(r(j), 1), vertices(r(j), 2));
            Vector3d n = Vector_Normal.row(r(j)).normalized();
            Vector4d normal (n(0), n(1), n(2), 0);
            v.normal = normal;
            vertex_attributes.push_back(v);
        }
    }

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

int main(int argc, char *argv[])
{
    setup_scene();

    int W = H * aspect_ratio;
    Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> frameBuffer(W, H);
    vector<uint8_t> image;

    simple_render(frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("simple.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);
    frameBuffer = Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic>(W, H);

    wireframe_render(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("wireframe.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);
    frameBuffer = Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic>(W, H);

    flat_shading(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("flat_shading.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);
    frameBuffer = Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic>(W, H);

    pv_shading(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("pv_shading.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);
    frameBuffer = Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic>(W, H);


    // add the animation
    int delay = 4;
    GifWriter g;
    GifBegin(&g, "wireframe.gif", frameBuffer.rows(), frameBuffer.cols(), delay);
    for (float i = 0; i < 2 * M_PI; i += 0.05)
    {
        frameBuffer.setConstant(FrameBufferAttributes());
        wireframe_render(i, frameBuffer);
        framebuffer_to_uint8(frameBuffer, image);
        GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
    }
    GifEnd(&g);

    GifBegin(&g, "flat_shading.gif", frameBuffer.rows(), frameBuffer.cols(), delay);
    for (float i = 0; i < 2 * M_PI; i += 0.05)
    {
        frameBuffer.setConstant(FrameBufferAttributes());
        flat_shading(i, frameBuffer);
        framebuffer_to_uint8(frameBuffer, image);
        GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
    }
    GifEnd(&g);

    GifBegin(&g, "pv_shading.gif", frameBuffer.rows(), frameBuffer.cols(), delay);
    for (float i = 0; i < 2 * M_PI; i += 0.05)
    {
        frameBuffer.setConstant(FrameBufferAttributes());
        pv_shading(i, frameBuffer);
        framebuffer_to_uint8(frameBuffer, image);
        GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
    }
    GifEnd(&g);

    return 0;
}