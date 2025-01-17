#include "window.h"

#include <png.h>
#include <glm/gtc/type_ptr.hpp>

#include "../util/physical_constants.h"

int openWindow(const char * vertex_file_path, const char * fragment_file_path) {
    int size[2] = {1000, 1000};
    return openWindowSized(vertex_file_path, fragment_file_path, size);
}

int openWindowSized(const char * vertex_file_path, const char * fragment_file_path, int size[2]) {
    // Initialise GLFW
    if( !glfwInit() )
    {
        fprintf( stderr, "Failed to initialize GLFW\n" );
        return -1;
    }

    glfwWindowHint(GLFW_SAMPLES, 0); // No antialiasing
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4); // We want OpenGL 4.5
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // We don't want the old OpenGL

    // Open a window and create its OpenGL context
    window = glfwCreateWindow( size[1], size[0], "Proton Radiography", NULL, NULL);
    if( window == NULL ){
        fprintf( stderr, "Failed to open window. This program requires OpenGL 4.5\n" );
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window); // Initialize GLEW
    glewExperimental=true; // Needed in core profile
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        return -1;
    }

    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    glfwSwapInterval(0);

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_PROGRAM_POINT_SIZE);

    programID = LoadShaders(vertex_file_path, fragment_file_path);

    MatrixID = glGetUniformLocation(programID, "MVP");

    return 0;
}

int setupMatrix(double max_beta, FieldStructure* field) {
    Projection = glm::ortho(-max_beta * c, max_beta * c, -max_beta * c, max_beta * c, -c, c);

    // Camera matrix
    View = glm::lookAt(
        glm::vec3(0, 0, 0), // Camera is at the origin of velocity space
        glm::vec3(0, 0, c), // and looks at +c
        glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
    );

    // Model matrix : an identity matrix (model will be at the origin)
    double matrix[16] = {
        field->xaxis.x, field->yaxis.x, field->zaxis.x, 0,
        field->xaxis.y, field->yaxis.y, field->zaxis.y, 0,
        field->xaxis.z, field->yaxis.z, field->zaxis.z, 0,
        0, 0, 0, 1
    };

    Model = glm::make_mat4(matrix);
    // Our ModelViewProjection : multiplication of our 3 matrices
    mvp = Projection * View * Model; // Remember, matrix multiplication is the other way around

    printf("[%e, %e, %e, %e;\n %e, %e, %e, %e;\n %e, %e, %e, %e;\n %e, %e, %e, %e]",
            mvp[0][0], mvp[0][1], mvp[0][2], mvp[0][3],
            mvp[1][0], mvp[1][1], mvp[1][2], mvp[1][3],
            mvp[2][0], mvp[2][1], mvp[2][2], mvp[2][3],
            mvp[3][0], mvp[3][1], mvp[3][2], mvp[3][3]
        );
    return 0;
}

int setupMatrix(ParticleSource* source, FieldStructure* field, ParticleDetector* detector) {
    double sourceDistance = source->distance;
    double openingAngle = source->divergence;
    double detectorDistance = detector->distance;

    Projection = glm::perspective(2 * openingAngle, 1.0, sourceDistance+field->min_z, sourceDistance + detectorDistance*1.01);

    // Camera matrix
    View = glm::lookAt(
        glm::vec3(0, 0, -sourceDistance), // Camera is at (0, 0, -d), in World Space
        glm::vec3(0, 0, detectorDistance), // and looks at the origin
        glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
    );

    // Model matrix : an identity matrix (model will be at the origin)
    double matrix[16] = {
        field->xaxis.x, field->yaxis.x, field->zaxis.x, 0,
        field->xaxis.y, field->yaxis.y, field->zaxis.y, 0,
        field->xaxis.z, field->yaxis.z, field->zaxis.z, 0,
        0, 0, 0, 1
    };

    Model = glm::make_mat4(matrix);
    // Our ModelViewProjection : multiplication of our 3 matrices
    mvp = Projection * View * Model; // Remember, matrix multiplication is the other way around

    return 0;
}

int setupBuffers(double pos[], bool running[], long N) {
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    // Generate 2 buffers, put the resulting identifier in vertexbuffer
    glGenBuffers(2, vertexbuffers);

    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffers[0]);
    glBufferData(GL_ARRAY_BUFFER, 3*N*sizeof(double), pos, GL_STREAM_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffers[1]);
    glBufferData(GL_ARRAY_BUFFER, N*sizeof(bool), running, GL_STREAM_DRAW);

    return 0;
}

int updateBuffers(double pos[], bool running[], long N) {
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffers[0]);
    glBufferSubData(GL_ARRAY_BUFFER, 0, 3*N*sizeof(double), pos);

    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffers[1]);
    glBufferSubData(GL_ARRAY_BUFFER, 0, N*sizeof(bool), running);

    return 0;
}

int setupBuffersTexRectangle(float* tbo_data, int pixels[2]) {
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    // Generate 2 buffers, put the resulting identifier in vertexbuffer
    glGenBuffers(2, vertexbuffers);

    float verts[12] = {-1, -1, 0, -1, 1, 0, 1, -1, 0, 1, 1, 0};
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffers[0]);
    glBufferData(GL_ARRAY_BUFFER, 3*4*sizeof(float), verts, GL_STATIC_DRAW);

    //glBindBuffer(GL_PIXEL_UNPACK_BUFFER, vertexbuffers[1]);
    //glBufferData(GL_PIXEL_UNPACK_BUFFER, n_pixels*sizeof(double), tbo_data, GL_STREAM_DRAW);
    glGenTextures(1, &TextureID);
    glBindTexture(GL_TEXTURE_2D, TextureID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, pixels[1], pixels[0], 0, GL_RED, GL_FLOAT, (GLvoid*)tbo_data);

    //glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);

    TexUniformID = glGetUniformLocation(programID, "tex");

    return 0;
}

int updateBuffersTexRectangle(float* tbo_data, int pixels[2]) {
    glBindTexture(GL_TEXTURE_2D, TextureID);
    //glBindBuffer(GL_PIXEL_UNPACK_BUFFER, vertexbuffers[1]);
    //glBufferSubData(GL_PIXEL_UNPACK_BUFFER, 0, n_pixels*sizeof(float), tbo_data);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, pixels[1], pixels[0], GL_RED, GL_FLOAT, tbo_data);

    //glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);

    return 0;
}

bool save_png_libpng(const char *filename, uint8_t *pixels, int w, int h) {
    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (!png)
        return false;

    png_infop info = png_create_info_struct(png);
    if (!info) {
        png_destroy_write_struct(&png, &info);
        return false;
    }

    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        png_destroy_write_struct(&png, &info);
        return false;
    }

    png_init_io(png, fp);
    png_set_IHDR(png, info, w, h, 8 /* depth */, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_colorp palette = (png_colorp)png_malloc(png, PNG_MAX_PALETTE_LENGTH * sizeof(png_color));
    if (!palette) {
        fclose(fp);
        png_destroy_write_struct(&png, &info);
        return false;
    }
    png_set_PLTE(png, info, palette, PNG_MAX_PALETTE_LENGTH);
    png_write_info(png, info);
    png_set_packing(png);

    png_bytepp rows = (png_bytepp)png_malloc(png, h * sizeof(png_bytep));
    for (int i = 0; i < h; ++i)
        rows[i] = (png_bytep)(pixels + (h - i) * w * 3);

    png_write_image(png, rows);
    png_write_end(png, info);
    png_free(png, palette);
    png_destroy_write_struct(&png, &info);

    fclose(fp);
    delete[] rows;
    return true;
}

int draw(int N, int name_length, char* name, bool wait) {
    do {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glUseProgram(programID);

        glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &mvp[0][0]);

        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, vertexbuffers[0]);
        glVertexAttribPointer(
            0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
            3,                  // size
            GL_DOUBLE,          // type
            GL_FALSE,
            0,                  // stride
            (void*)0            // array buffer offset
        );

        glEnableVertexAttribArray(1);
        glBindBuffer(GL_ARRAY_BUFFER, vertexbuffers[1]);
        glVertexAttribIPointer(
            1,                  // attribute 1. No particular reason for 1, but must match the layout in the shader.
            1,                  // size
            GL_UNSIGNED_BYTE,   // type
            0,                  // stride
            (void*)0            // array buffer offset
        );

        // Draw the triangle !
        glDrawArrays(GL_POINTS, 0, N);
        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);

        // Swap buffers
        glfwSwapBuffers(window);
        glfwPollEvents();
    }   while( wait && glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
               glfwWindowShouldClose(window) == 0 );
    if ( name_length > 0 ) {
        int w = 1000;
        int h = 1000;
        uint8_t *pixels = new uint8_t[w * h * 3];
        // copy pixels from screen
        glMemoryBarrier(GL_ALL_BARRIER_BITS);
        glReadBuffer(GL_FRONT);
        glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid *)pixels);

        // save the image
        int err = save_png_libpng(name, pixels, w, h);
        if (err)
           printf("Done\n");
        else
           printf("Failed\n");
    }
    return 0;
}

int drawTexRectangle() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glUseProgram(programID);


    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffers[0]);
    glVertexAttribPointer(
        0,          // attribute 0. No particular reason for 0, but must match the layout in the shader.
        3,          // size
        GL_FLOAT,   // type
        GL_FALSE,
        0,          // stride
        (void*)0    // array buffer offset
    );

    glActiveTexture(GL_TEXTURE0);
    glUniform1i(TexUniformID, 0);

    // Draw the quad !
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    glDisableVertexAttribArray(0);

    // Swap buffers
    glfwSwapBuffers(window);
    glfwPollEvents();

    return 0;
}


