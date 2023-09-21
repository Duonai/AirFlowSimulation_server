#define _CRT_SECURE_NO_WARNINGS

#include <cstring>
#include <iostream>
#include "Packet.h"

namespace tcp_boost
{
    int Packet::_buffer_length = Packet::default_buffer_length;


    int
    Packet::get_buffer_length()
    {
        return _buffer_length;
    }

    void
    Packet::set_buffer_length(int value)
    {
        _buffer_length = value;
    }

    Packet::Packet(size_t size)
        : m_position(header_length),
          m_body_length(0)
    {
        _buffer_length = size;

        m_buffer = new char[_buffer_length]();
    }

    Packet::Packet(const unsigned char *buffer, size_t size)
        : m_position(header_length),
          m_body_length(0)
    {
        _buffer_length = size;
            
        // deep copy
        m_buffer = (char*)std::memcpy(new char[_buffer_length], buffer, _buffer_length);
        decode_body_length();
    }

    Packet::Packet(const Packet &orig)
        : m_position(orig.get_position()),
          m_body_length(orig.get_body_length())
    {
        // deep copy
        m_buffer = (char*)std::memcpy(new char[_buffer_length], orig.get_buffer(), _buffer_length);
    }

    Packet::~Packet()
    {
        if (m_buffer)
            delete[] m_buffer;
    }

    void
    Packet::record_body_length()
    {
        char *header = (char*)&m_body_length;

        m_buffer[0] = m_body_length >> 24 & 0xFF;
        m_buffer[1] = m_body_length >> 16 & 0x00FF; //here
        m_buffer[2] = m_body_length >> 8 & 0x0000FF;
        m_buffer[3] = m_body_length & 0x000000FF;
        

        //std::memcpy(m_buffer, header, header_length);
    }

    void
    Packet::decode_body_length()
    {
        unsigned char header[header_length];
        std::memcpy(header, m_buffer, header_length);

        m_body_length = (header[0] << 24) | (header[1] << 16) | (header[2] << 8) | header[3];
    }

    char*
    Packet::get_buffer() const
    {
        return m_buffer;
    }

    char*
    Packet::get_body() const
    {
        return m_buffer + header_length;
    }

    int
    Packet::get_position() const
    {
        return m_position;
    }

    int
    Packet::get_total_length() const
    {
        return m_body_length + header_length;
    }

    int
    Packet::get_body_length() const
    {
        return m_body_length;
    }

    char
    Packet::pop_byte()
    {
        size_t size = 1;
        char *data = new char[size];
        read_buffer(data, size);
        return *data;
    }

    bool
    Packet::pop_bool()
    {
        char data = pop_byte();
        return data != 0;
    }

    int16_t
    Packet::pop_int16()
    {
        size_t size = sizeof(int16_t);
        char *data = new char[size];
        read_buffer(data, size);
        return *((int16_t*)data);
    }

    uint16_t
        Packet::pop_uint16()
    {
        size_t size = sizeof(uint16_t);
        char* data = new char[size];
        read_buffer(data, size);
        return *((uint16_t*)data);
    }

    int32_t
    Packet::pop_int32()
    {
        size_t size = sizeof(int32_t);
        char *data = new char[size];
        read_buffer(data, size);
        return *((int32_t*)data);
    }

    int64_t
    Packet::pop_int64()
    {
        size_t size = sizeof(int64_t);
        char *data = new char[size];
        read_buffer(data, size);
        return *((int64_t*)data);
    }

    int32_t
        Packet::pop_uint32()
    {
        size_t size = sizeof(uint32_t);
        char* data = new char[size];
        read_buffer(data, size);
        return *((uint32_t*)data);
    }

    float
    Packet::pop_single()
    {
        size_t size = sizeof(float);
        char *data = new char[size];
        read_buffer(data, size);
        unsigned long value;

        value = ((data[3] & 0xFF) << 24);
        value += ((data[2] & 0xFF) << 16);
        value += ((data[1] & 0xFF) << 8);
        value += data[0] & 0xFF;

        return *((float*)&value);
        //float value = 0;
        //memcpy(&value, &data, sizeof(float));
        //return ;// value; //float(data[0] << 24 | data[1] << 16 | data[2] << 8 | data[3]);
    }

    double
    Packet::pop_double()
    {
        size_t size = sizeof(double);
        char *data = new char[size];
        read_buffer(data, size);
        return *((double*)data);
    }

    char*
    Packet::pop_byte_array()
    {
        size_t size = pop_int64();
        char *data = new char[size];
        read_buffer(data, size);
        return data;
    }

    std::string
    Packet::pop_string()
    {
        char *data = pop_byte_array();
        return std::string(data);
    }

    void
    Packet::push_byte(char data)
    {
        size_t size = 1;
        char *temp = &data;
        write_buffer(temp, size);
    }

    void
    Packet::push_bool(bool data)
    {
        push_byte((char)data);
    }

    void
    Packet::push_int16(int16_t data)
    {
        size_t size = sizeof(int16_t);
        char *temp = (char*)&data;
        write_buffer(temp, size);
    }

    void
        Packet::push_uint16(uint16_t data)
    {
        size_t size = sizeof(uint16_t);
        char* temp = (char*)&data;
        write_buffer(temp, size);
    }

    void
    Packet::push_int32(int32_t data)
    {
        size_t size = sizeof(int32_t);
        char *temp = (char*)&data;
        write_buffer(temp, size);
    }

    void
    Packet::push_int64(int64_t data)
    {
        size_t size = sizeof(int64_t);
        char *temp = (char*)&data;
        write_buffer(temp, size);
    }

    void
    Packet::push_single(float data)
    {
        size_t size = sizeof(float);
        char *temp = (char*)&data;
        write_buffer(temp, size);
    }

    void
    Packet::push_double(double data)
    {
        size_t size = sizeof(double);
        char *temp = (char*)&data;
        write_buffer(temp, size);
    }

    void
    Packet::push_byte_array(char *data, size_t length)
    {
        //push_int64(length);
        write_buffer(data, length);
    }

    void
    Packet::push_string(std::string data)
    {
        size_t size = data.length() + 1;
        char *temp = strcpy(new char[size], data.c_str());
        push_byte_array(temp, size);
    }

    char*
    Packet::read_buffer(char *dest, size_t length)
    {
        // overflow
        if (m_position + length > get_total_length())
        {
            std::cout << "Packet reading warning : You are trying to read buffer over the total length.\n";
        }
        else if (m_position + length > _buffer_length)
        {
            std::cout << "Packet reading failed : Tryed to read buffer over the buffer length.\n";
            return nullptr;
        }

        std::memcpy(dest, m_buffer + m_position, length);
        m_position += length;
       
        return dest;
    }

    void
    Packet::write_buffer(char *src, size_t length)
    {
        if (m_position + length > _buffer_length)
        {
            std::cout << "Packet writing failed : Tryed to write buffer over the buffer length.\n";
            return;
        }

        std::memcpy(m_buffer + m_position, src, length);
        m_position += length;
        m_body_length += length;
        record_body_length();
    }
}