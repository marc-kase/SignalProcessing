package src;


import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class ListStreamTest {
    public static void main(String[] args) {
        String[] arr = {"1","2","3","4","5","6","7","8","9"};
        List<String> list = Arrays.asList(arr);
        List res = list.stream().limit(3).collect(Collectors.toList());
        res.stream().forEach(System.out::println);
    }
}
