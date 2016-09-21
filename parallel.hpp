#pragma once
#include <functional>
#include <vector>
#include <type_traits>
#include <mpi.h>

// TODO: implement async

class Parallel {
    const int kRoot = 0;
    int rank_;
    int size_;
    
public:
    Parallel();

    int rank();

    int size();
    
    bool is_root();
    // Executes something on root rank
    void ExecRoot(std::function<void()> target);
    
    template<typename T> auto Reduce(T value, std::function<T(std::vector<T>&)> reduce) -> 
    std::enable_if_t<std::is_trivially_copyable<T>::value, T>;
    template<typename T> auto ReduceToAll(T value, std::function<T(std::vector<T>&)> reduce) -> 
    std::enable_if_t<std::is_trivially_copyable<T>::value, T>;
    template<typename T> auto ReduceRootToAll(T value, std::function<T(std::vector<T>&)> reduce) -> 
    std::enable_if_t<std::is_trivially_copyable<T>::value, T>;

    template<typename T> auto Reduce(T value, std::function<T(std::vector<T>&)> reduce) -> 
    std::enable_if_t<std::is_trivially_copyable<typename T::value_type>::value, std::enable_if_t<std::is_same<std::vector<typename T::value_type>, T>::value, T>>;
    template<typename T> auto ReduceToAll(T value, std::function<T(std::vector<T>&> reduce) -> 
    std::enable_if_t<std::is_trivially_copyable<typename T::value_type>::value, std::enable_if_t<std::is_same<std::vector<typename T::value_type>, T>::value, T>>;


    template<typename T> auto Gather(T value) ->
    std::enable_if_t<std::is_trivially_copyable<T>::value, std::vector<T>>;

    template<typename T> auto AllGather(T value) ->
    std::enable_if_t<std::is_trivially_copyable<T>::value, std::vector<T>>;

    void Barrier();
};

template<typename T> auto Parallel::Reduce(T value, std::function<T(std::vector<T>&)> reduce) -> 
std::enable_if_t<std::is_trivially_copyable<T>::value, T> {
    std::vector<T> data_buffer;
    if(is_root()) {
        data_buffer.resize(size());
    }

    MPI_Gather(&value, sizeof(T), MPI_BYTE, data_buffer.data(), sizeof(T), MPI_BYTE, kRoot, MPI_COMM_WORLD);

    if(is_root()) {
        return reduce(data_buffer);
    }else {
        return {};
    }
}

template<typename T> auto Parallel::ReduceRootToAll(T value, std::function<T(std::vector<T>&)> reduce) -> 
std::enable_if_t<std::is_trivially_copyable<T>::value, T> {
    T reducedvalue = Reduce(value, reduce);

    MPI_BCAST(&reducedvalue, sizeof(T), MPI_BYTE, kRoot, MPI_COMM_WORLD);
    return reducedvalue;
}

template<typename T> auto Parallel::ReduceToAll(T value, std::function<T(std::vector<T>&)> reduce) -> 
std::enable_if_t<std::is_trivially_copyable<T>::value, T> {
    std::vector<T> data_buffer(sizeof(T) * size());

    MPI_Allgather(&value, sizeof(T), MPI_BYTE, data_buffer.data(), sizeof(T), MPI_BYTE, MPI_COMM_WORLD);

    return reduce(data_buffer);
}

template<typename T> auto Parallel::Reduce(T value, std::function<T(std::vector<T>&)> reduce) -> 
std::enable_if_t<std::is_trivially_copyable<typename T::value_type>::value, std::enable_if_t<std::is_same<std::vector<typename T::value_type>, T>::value, T>> {
    
    if(is_root()) {
        MPI_Status status;
        MPI_Message message;
        int message_size;
        std::vector<T> data_buffers(size());
        data_buffers.front() = value;
        
        // TODO: make this async
        for(std::size_t k = 1; k < size(); ++k) {
            MPI_Mprobe(k, 0, MPI_COMM_WORLD, &message, &status);
            MPI_Get_count(&status, MPI_BYTE, &message_size);
            data_buffers[k].resize(message_size);
            MPI_Mrecv(data_buffers[k].data(), message_size, MPI_BYTE, &message, &status);
        }

        return reduce(data_buffers);
    }else {
        MPI_Bsend(value.data(), value.size() * sizeof(value.value_type), MPI_BYTE, kRoot, 0, MPI_COMM_WORLD);
        return { };
    }
}