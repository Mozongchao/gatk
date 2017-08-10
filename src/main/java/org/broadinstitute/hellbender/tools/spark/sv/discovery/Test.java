package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import java.util.Spliterator;
import java.util.function.Consumer;

/**
 * Created by valentin on 8/9/17.
 */
public class Test {
    public static void main(final String[] args) {



    }

    class Sp<T> implements Spliterator<T> {

        @Override
        public boolean tryAdvance(Consumer<? super T> action) {
            return false;
        }

        @Override
        public Spliterator<T> trySplit() {
            return null;
        }

        @Override
        public long estimateSize() {
            return 100;
        }

        @Override
        public int characteristics() {
            return 0;
        }
    }
}
