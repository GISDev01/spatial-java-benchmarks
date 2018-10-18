package com.gisdev01;

import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import static java.util.stream.Collectors.*;


public class Main {

    public static void main(String[] args) throws Exception {

        Map<String, Integer> mapStrIntTest = new HashMap<>();
        mapStrIntTest.put("KeyA", 120);
        mapStrIntTest.put("KeyB", 150);
        mapStrIntTest.put("KeySuperHigh", 1150);
        mapStrIntTest.put("KeySuperLow", 10);
        mapStrIntTest.put("KeyC", 100);

        System.out.println("Map pre-sort: " + mapStrIntTest);

        // Sort the map in decreasing order of value
        Map<String, Integer> sorterDescMap = mapStrIntTest
                .entrySet()
                .stream()
                .sorted(Collections.reverseOrder(Map.Entry.comparingByValue()))
                .collect(
                        toMap(Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e2,
                                LinkedHashMap::new));

        System.out.println("map after sorting by values in descending order: "
                + sorterDescMap);

        for (Map.Entry<String, Integer> entry : sorterDescMap.entrySet()) {
            String key = entry.getKey();
            Integer value = entry.getValue();
            System.out.println(key + ": " + value);
        }
    }

}
