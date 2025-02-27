
#include "vmec_output.h"
#include "vmec_common.h"

// External Includes
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"
#pragma clang diagnostic pop

#include <ImfOutputFile.h>
#include <ImfHeader.h>
#include <ImfArray.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>

namespace fsys = std::filesystem;

void image_generate(pix_realNumber* pixels_in, pix_RGB* pixels_out)
{
    realNumber color_temp;
    for(int y=0;y<settings.image_height;++y)
    {
        for(int x=0;x<settings.image_width;++x)
        {
            color_temp=std::min(std::max(float(pixels_in[x+y*settings.image_width].R),0.0f),1.0f);
            pixels_out[x+y*settings.image_width].R=static_cast<unsigned char>(255.0f*color_temp);
            color_temp=std::min(std::max(float(pixels_in[x+y*settings.image_width].G),0.0f),1.0f);
            pixels_out[x+y*settings.image_width].G=static_cast<unsigned char>(255.0f*color_temp);
            color_temp=std::min(std::max(float(pixels_in[x+y*settings.image_width].B),0.0f),1.0f);
            pixels_out[x+y*settings.image_width].B=static_cast<unsigned char>(255.0f*color_temp);
        }
    }
}


void write_binary_RBG(pix_RGB* pixels_all, int frame_number) {
    char filename[30];
    std::memset(filename, 0, 30);
    std::sprintf(filename, "Render%03d.bin", frame_number);

    fpath file_path = fpath(settings.render_binary_dir) / filename;

    int img_width = settings.image_width;
    int img_height = settings.image_height;
    std::ofstream fs(filename, std::ios::out | std::ios::binary | std::ios::trunc);
    fs.write(reinterpret_cast<const char*>(&img_width), sizeof(int));
    fs.write(reinterpret_cast<const char*>(&img_height), sizeof(int));
    fs.write(reinterpret_cast<const char*>(pixels_all), settings.image_height*settings.image_width*3);
    fs.close();
    vmec_log("Wrote binary to file: ", file_path);
}



void write_binary_realNumber(pix_realNumber* pixels_all, int frame_number) {
    char filename[30];
    std::memset(filename, 0, 30);
    std::sprintf(filename, "Render%03d.bin", frame_number);

    fpath file_path = fpath(settings.render_binary_dir)  / filename;

    int img_width = settings.image_width;
    int img_height = settings.image_height;
    std::ofstream fs(file_path, std::ios::out | std::ios::binary | std::ios::trunc);
    fs.write(reinterpret_cast<const char*>(&img_width), sizeof(int));
    fs.write(reinterpret_cast<const char*>(&img_height), sizeof(int));
    fs.write(reinterpret_cast<const char*>(pixels_all), settings.image_height*settings.image_width*3*sizeof(realNumber));
    fs.close();
    vmec_log("Wrote vmec_binary to file: ", file_path);
}

void write_hdr_realNumber(pix_realNumber* pixels_all, int frame_number)
{
    char filename[30];
    std::memset(filename, 0, 30);
    std::sprintf(filename, "Render%03d.hdr", frame_number);

    fpath file_path = fpath(settings.render_png_dir)  / filename;

    i32 width = settings.image_width;
    i32 height = settings.image_height;
    std::ofstream fs(file_path, std::ios::out | std::ios::binary | std::ios::trunc);
    pixel_f32* buffer = allocate<pixel_f32>(width * height * 3);
    static_assert( sizeof(realNumber) == 8);

     for (int i=0; i< width * height; ++i)
     {
         buffer[i].R = (f32)pixels_all[i].R;
         buffer[i].G = (f32)pixels_all[i].G;
         buffer[i].B = (f32)pixels_all[i].B;
     }
     std::string msvc_copy = file_path.string();
     bool success = stbi_write_hdr(msvc_copy.c_str(), width, height, 3, (f32*)buffer);
     if (success )
     { vmec_log("Wrote HDR to file: ", file_path); }
     // Cleanup
     free(buffer);
}

void
write_exr_realNumber(pix_realNumber* pixels_all, int frame_number)
{
    char filename[30];
    std::memset(filename, 0, 30);
    std::sprintf(filename, "Render%03d.exr", frame_number);

    fpath file_path = fpath(settings.render_exr_dir)  / filename;

    i32 width = settings.image_width;
    i32 height = settings.image_height;

    pixel_f32* buffer = allocate<pixel_f32>(height * width * 3);

     for (int i=0; i< height * width; ++i)
     {
         buffer[i].R = (f32)pixels_all[i].R;
         buffer[i].G = (f32)pixels_all[i].G;
         buffer[i].B = (f32)pixels_all[i].B;
     }

    Imf::Header exr_header(width, height);
    Imf::ChannelList& channels = exr_header.channels();

    channels.insert( "R", Imf::Channel(Imf::FLOAT) );
    channels.insert( "G", Imf::Channel(Imf::FLOAT) );
    channels.insert( "B", Imf::Channel(Imf::FLOAT) );

    std::string msvc_copy = file_path.string();
    Imf::OutputFile exr_file(msvc_copy.c_str(), exr_header);
    Imf::FrameBuffer framebuffer;

    i32 pixel_stride = (sizeof(pixel_f32));
    i32 row_stride = (pixel_stride * width);

    f32* buffer_channels = (f32*)buffer;
    auto red_channel = Imf::Slice( Imf::FLOAT, (char*)(buffer_channels+ 0), pixel_stride, row_stride);
    auto green_channel = Imf::Slice( Imf::FLOAT, (char*)(buffer_channels+ 1), pixel_stride, row_stride);
    auto blue_channel = Imf::Slice( Imf::FLOAT, (char*)(buffer_channels+ 2), pixel_stride, row_stride);
    // auto alpha_channel = Imf::Slice( Imf::DOUBLE, (u8*)buffer + 3, pixel_stride, row_stride);

    framebuffer.insert( "R", red_channel );
    framebuffer.insert( "G", green_channel );
    framebuffer.insert( "B", blue_channel );
    // framebuffer.insert( "A", alpha_channel );

    exr_file.setFrameBuffer( framebuffer );
    exr_file.writePixels( height );
    vmec_log("Wrote EXR to file: ", file_path);

    // Cleanup
    free(buffer);
}
