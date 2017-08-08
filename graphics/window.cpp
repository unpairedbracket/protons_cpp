#include "window.h"

int openWindow(const char * vertex_file_path,const char * fragment_file_path) {
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
    window = glfwCreateWindow( 1000, 1000, "Particle Tracer", NULL, NULL);
    if( window == NULL ){
        fprintf( stderr, "Failed to open window. This program requires OpenGL 4.\n" );
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

    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glEnable(GL_PROGRAM_POINT_SIZE);
    
    programID = LoadShaders(vertex_file_path, fragment_file_path);

    MatrixID = glGetUniformLocation(programID, "MVP");

    return 0;
}

int setupMatrix(double sourceDistance, double detectorDistance, double openingAngle) {
	// Projection matrix : 45° Field of View, 4:3 ratio, display range : 0.0 unit <-> 1 unit
    std::cout << "Angle: " << openingAngle << std::endl;
	Projection = glm::perspective(2 * openingAngle, 1.0, sourceDistance, sourceDistance + detectorDistance);

	std::cout << Projection[0][0] << "," << Projection[1][0] << "," << Projection[2][0] << "," << Projection[3][0] << std::endl;
	std::cout << Projection[0][1] << "," << Projection[1][1] << "," << Projection[2][1] << "," << Projection[3][1] << std::endl;
	std::cout << Projection[0][2] << "," << Projection[1][2] << "," << Projection[2][2] << "," << Projection[3][2] << std::endl;
	std::cout << Projection[0][3] << "," << Projection[1][3] << "," << Projection[2][3] << "," << Projection[3][3] << std::endl;

	// Camera matrix
	View = glm::lookAt(
		glm::vec3(0, 0, -sourceDistance), // Camera is at (0, 0, -d), in World Space
		glm::vec3(0, 0, detectorDistance), // and looks at the origin
		glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
	);

	// Model matrix : an identity matrix (model will be at the origin)
	Model = glm::mat4(1.0f);
	// Our ModelViewProjection : multiplication of our 3 matrices
	mvp = Projection * View * Model; // Remember, matrix multiplication is the other way around
	
	std::cout << mvp[0][0] << "," << mvp[1][0] << "," << mvp[2][0] << "," << mvp[3][0] << std::endl;
	std::cout << mvp[0][1] << "," << mvp[1][1] << "," << mvp[2][1] << "," << mvp[3][1] << std::endl;
	std::cout << mvp[0][2] << "," << mvp[1][2] << "," << mvp[2][2] << "," << mvp[3][2] << std::endl;
	std::cout << mvp[0][3] << "," << mvp[1][3] << "," << mvp[2][3] << "," << mvp[3][3] << std::endl;

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
    
int draw(int N, bool wait) {
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
    
    return 0;
}


